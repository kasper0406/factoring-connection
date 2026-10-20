# heurnfs: the heuristic L[1/3] frontier and the project's devices

Agent `heurnfs`, phase 14, 2026-10-20.
- Code: `src/heurnfs/models.py` (every cost program and closed form below), `src/heurnfs/optim.py` (seeded
  Nelder–Mead; floating point is used only to minimise constant programs, never in an operational routine).
- Tests: `tests/test_heurnfs.py` (13 tests: closed forms, numeric optima against them, the exact height LP,
  the exponent barriers).
- GPT-6.1-sol: blind round `codex/heurnfs_1.{md,out}`. It reached the same verdict independently (§6).
- Reviews: one independent Opus critic (model and proofs) and one Opus literature critic (§7).
- No experiment was run, so there is no pre-registration. The optional smoothness-yield experiment was dropped
  because no variant survived the analysis and so there was nothing to test.

## Verdict

**No heuristic improvement was found.** Nothing found here beats Coppersmith's MNFS constant
c_M = ((92+26√13)/27)^{1/3} = 1.9019 for one general N. Nothing gives a heuristic exponent below 1/3.

- **Every device leaves the constant unchanged or makes it worse.** These are Dixon descent, Theorem K
  splits, recursion, the N·w multiplier, the 2D-random rational side and batch smoothness.
  - The devices were built to make *proofs* work: fruitfulness, relation supply, shape tests. Classical
    heuristic H already assumes all of that away.
  - None of them touches the three things that set the constant:
    1. the norm-size ledger 2x/d + d·log E;
    2. the sieve cost E²;
    3. linear algebra (LA) at ω = 2.
- **Why the MNFS optimum is where it is (model-PROVED, i.e. exact within heuristic H).**
  - Three constraints bind: sieve = LA = 1.9019, together with the relation count.
  - The algebraic-test cost has slack: 1.768 against 1.902. So cheaper smoothness tests (ECM, batch,
    recursive NFS) cannot help.
  - Block-angular LA cannot help either, even optimistically: it still gives 1.9019.
- **Precise "what it would take" (HYPOTHETICAL).**
  - (a) Enumerate the rational-smooth pairs of the box at output-sensitive cost, i.e. without paying E².
    MNFS would then fall to 1.6386, the factory's per-N constant.
  - (b) Sparse LA with exponent ω ≤ 1.5. MNFS would then fall to 1.8597. GNFS stays at 1.9230, because
    its sieve alone binds.
  - Neither (a) nor (b) is known, and no project device supplies either.
- **The exponent barrier (model-PROVED).**
  - Under the counting heuristic G for polynomial heights, every two-polynomial selection obeys
    log X ≥ 2x/D (1+o(1)) + D·log E. This holds for any degrees (d₁, d₂), any sieving dimension t, and tower
    (TNFS-type) rings.
  - So α ≥ 1/3. For general N, the 2x/D term is the counting barrier. SNFS beats it only because its N
    has special form (κ = 1, constant 1.526).
  - The two r20 long shots (toric entropy–height separation; Witt carries) each need a missing arithmetic
    construction (§5). Neither has a constant that can be computed.

Status labels:
- **KNOWN** means published, with a citation.
- **model-PROVED** means an exact consequence of the stated heuristic cost model, checked in code.
- **HEURISTIC** means it rests on an assumption that is named in place.
- **HYPOTHETICAL** means it is conditional on an algorithm that does not exist.

---

## 1. The heuristic model and the classical optimisations

**Heuristic H [assumption].**
- Let x = log N, and write L[c] for L_N[1/3, c].
- Polynomial values of size exp(s·x^{2/3}(log x)^{1/3}) are B-smooth, with B = L[β], with probability
  L[−s/(3β)]. This is Canfield–Erdős–Pomerance applied to polynomial values.
- Joint smoothness events multiply, up to L^{o(1)}.
- One smoothness test costs L^{o(1)} (ECM: L_B[1/2, √2] = L^{o(1)}; or batch).
- Sparse LA on C columns of weight L^{o(1)} costs C^{2+o(1)} (block Wiedemann).
- Dependencies are fruitful, and square roots are cheap.

**Heuristic G [assumption].** For general N, a polynomial system with heights H_i and a prescribed root
mod N exists only if the lattice count (∏ H_i^{#coeffs}) / N^{#constraints} is ≥ 1 (Gaussian heuristic). SNFS
inputs violate G, and that is exactly why SNFS is faster.

**GNFS [KNOWN; derivation reproduced].**
- Parameters: d = δ(x/log x)^{1/3}, base-m with m ≈ N^{1/(d+1)}, box |a|, |b| ≤ E = L[e], B = L[β].
- Rational norm: E·N^{1/d}, size coefficient 1/δ.
- Algebraic norm: E^d·N^{1/d}, size δe + 1/δ.
- Relations must exceed the factor base: 2e − (δe + 2/δ)/(3β) ≥ β.
- Cost: max(2e [sieve], 2β [LA]).
- At e = β the constraint reads 3β² ≥ δβ + 2/δ. Minimising the right side at δ = √(2/β) gives β³ = 8/9.
- So **c = 2β = (64/9)^{1/3} = 1.9230**, with δ = 3^{1/3}.
- More generally, with κ/δ in place of 2/δ: **c(κ) = (32κ/9)^{1/3}** (`gnfs_closed`). κ = 1 gives SNFS 1.5263.

**MNFS [KNOWN, Coppersmith 1993; reproduced].**
- One rational side g. V = L[v] algebraic sides f + (c₁x + c₀)g. The rational factor base is L[β]; each
  algebraic factor base is L[β_a].
- The columns number s = max(β, β_a + v).
- Relations: 2e + v − 1/(3βδ) − (δe + 1/δ)/(3β_a) ≥ s.
- Cost: max(2e, 2e − 1/(3βδ) + v [algebraic tests], 2s).
- Eliminating variables in closed form gives the binding condition 2β − 1/(3βδ) ≥ 2√((δβ + 1/δ)/3), which
  yields **c_M = 1.9019**. The optimum is:

| e = β | δ | v | β_a | algebraic test cost |
|---|---|---|---|---|
| 0.95094 | 1.40175 | 0.12503 | 0.82591 | 1.7768 |

  - (`mnfs`. The numeric minimiser lands on a flat ridge, v ∈ [0.116, 0.25], all at cost 1.9019–1.9020.)
  - The rational penalty 1/(3βδ) = 0.2501 equals 2v (also in codex heurnfs_1).
  - V·B_a = B: the rational and algebraic column counts balance.

**Factorization factory [KNOWN, Coppersmith 1993].**
- Per-N cost max(2e − r, 2max(β, β_a)), where the rational table is shared and costs 2e.
- Our program: **1.6386 per N** after a shared precomputation of ≈ 2.006 (literature 2.007; our minimiser's
  ridge gives 2.0055–2.008).
- For a single fresh N the precomputation is chargeable. The cost is then ≥ 2.006 > c_M.

**Batch NFS [KNOWN, Bernstein–Lange 2014, ePrint 2014/921].**
- L^{1.704…} per integer, amortized over a batch of L[1/3, 0.5] integers.
- The measure is *circuit area–time*, not RAM time.
- It is not a single-N RAM constant, so it does not compete with c_M (codex heurnfs_1 §1).

**Other published variants, for completeness [KNOWN].**
- exTNFS (Kim–Barbulescu, 1.7472), MexTNFS (1.7114), SexTNFS (1.5263) and Sarkar–Singh's algorithms are
  DL algorithms in F_{p^n} with composite n.
- They exploit a known characteristic p and a subfield structure: the polynomial constraint is modulo
  p, not modulo Q = p^n. §4 Barrier B3 shows that for factoring the tower degree cancels.
- Lee–Venkatesan 2018 (arXiv 1805.08873) prove rigorously that their randomized NFS (and its MNFS form)
  produces congruences of squares at the GNFS (respectively MNFS) constant. Extracting a factor stays
  conditional on their Conjecture 7.1.
- Literature search for 2018–2026, by codex heurnfs_1 and our literature critic: no validated
  factoring constant below c_M. Schnorr's 2021 claims are unsubstantiated. Kirchner–Fouque 2020 is about
  parallel and area–time cost.

## 2. Table 1: constants per variant (all L_N[1/3, ·], single N, RAM time)

Every number is produced by `src/heurnfs/models.py::table()` (run: `.venv/bin/python -m src.heurnfs.models`).

| # | Variant | Constant | Status | Why |
|---|---|---|---|---|
| 1 | GNFS | **1.9230** | KNOWN | §1 |
| 2 | MNFS (Coppersmith) | **1.9019** | KNOWN | §1 |
| 3 | Factory, per N / precomp | 1.6386 / 2.006 | KNOWN | precomp chargeable for one N |
| 4 | SNFS (κ = 1), special N only | 1.5263 | KNOWN | violates heuristic G; not general N |
| 5 | Sieving dimension t = 3, 4 (GNFS) | 2.1165, 2.2013 | KNOWN (folklore); closed form (128(t−1)/(9t))^{1/3} | ledger 2(t−1)x/D + D·C/t; t = 2 optimal |
| 6 | Two degree-d polys (Montgomery/JL-type) | 1.9230 | model-PROVED | ledger x/d + 2d log E: same as GNFS after δ ↦ 2δ |
| 7 | Tower ring Z[ι], any η (TNFS-type) | = GNFS | model-PROVED | η cancels in the height LP (`tower_ledger`) |
| 8 | "Two-of-V" symmetric MNFS | 2.3962 | model-PROVED | both sides pay N^{1/(d+1)}·E^d |
| 9 | Multiple rational sides | not available / ≥ GNFS | PROVED (elementary) | §3.5 |
| 10 | N·w multiplier, w ≈ N^{1/d} in a₀ (ldbound's literal design) | GNFS 2.2013 = (32/3)^{1/3}; MNFS 2.1909 | model-PROVED | κ: 2 → 3 |
| 11 | N·w re-expanded in base (Nw)^{1/(d+1)}, w = N^{ω} | 1.9230·(1+ω)^{1/3} (GNFS) | model-PROVED | ω → 0 at growing degree: o(1) |
| 12 | 2D-random rational side g₁x − g₀ (pcest) | 1.9230 / 1.9019 | model-PROVED | identical norm sizes |
| 13 | Batch smoothness / ECM / recursive-NFS cofactorization | 1.9230 / 1.9019 | model-PROVED | test cost already L^{o(1)} and non-binding |
| 14 | Bounded large primes ≤ B^{O(1)} | 1.9230 / 1.9019 | KNOWN; model-PROVED | gain ρ(u−O(1))/ρ(u) = u^{O(1)} = L^{o(1)} |
| 15 | Descent rows, bottom bound (any descent tree) | ≥ 1.9230 (single field), ≥ c_M with MNFS bottom | model-PROVED | §3.1 |
| 16 | Dixon-descent rows, plain initial split (3^{1/3}) | ≥ 2.1144 | model-PROVED (lower bound) | §3.1 |
| 17 | Theorem K split + descent of L[2/3] primes in the main field | no finite L[1/3] constant | model-PROVED | §3.2 |
| 18 | MNFS + block-angular LA (optimistic) | 1.9019 | HYPOTHETICAL; model-PROVED no gain | LA structure does not bind |
| 19 | MNFS + output-sensitive rational-smooth enumeration | **1.6386** | HYPOTHETICAL | §4.4 (a) |
| 20 | Sparse LA exponent ω = 1, 1.5 / 2 / 2.373 | MNFS 1.8597 / 1.9019 / 1.9555; GNFS 1.9230 / 1.9230 / 1.9553 | HYPOTHETICAL (ω < 2) | §4.4 (b) |

## 3. The devices, one by one (Step 2)

### 3.1 Dixon descent / descent-style relation generation

**The device.** lemmaM Theorem DD:
- expose R = r² mod N;
- split R ≡ u/v by rational reconstruction;
- descend the primes of uv by special-q relations.

Rows built this way are fruitful by fibre symmetry, with probability 1 − 2^{1−ω}. Under heuristic H,
fruitfulness is *assumed*, so this benefit is worth nothing.

**Bottom bound [model-PROVED under H and assumption HD].**
- *HD (codex):* the estimates of H also hold for descent-selected children, with branching overhead
  L^{o(1)}.
- *Columns.* In a system of descent trees, every large prime is created once, in its parent, and
  removed once, by its child. So (#relations) − (#large primes) = #rows, and a dependency needs
  #rows > |FB| = L[β].
- *Each row.* Every row contains at least one bottom relation: a special-q relation with q ~ B whose
  other primes are all ≤ B.
- *Bottom relation cost.* Its q-lattice region must contain L[P] points, so 2e ≥ β + P with
  P = (δe + 2/δ)/(3β). Dividing by q = L[1/3]-size leaves the L[2/3]-scale norm coefficients unchanged.
- *Total.* ≥ L[β + P], and LA costs L[2β]. This is *exactly* the GNFS system, so the bound is 1.9230
  (`descent_rows()`). With MNFS-shared bottom relations it is c_M.
- *Conclusion.* Descent can never beat sieving.
  - Special-q lattice sieving is descent with tree depth 1.
  - Its equivalence to box sieving, up to o(1), is folklore.
- *Assumption.* Distinct rows use distinct bottom relations. If rows shared bottom relations, the
  column count would rise and the identity above would fail. We did not find a way to exploit sharing.

**Random-square rows [model-PROVED lower bound].**
- The plain initial split needs u, v ~ N^{1/2} to be L[2/3, y]-smooth, recognised by ECM. It costs
  1/(3y) + √(4y/3), minimised at y = 3^{−1/3} to **3^{1/3} = 1.442** (`init_split_closed`; this is
  Commeine–Semaev's first step).
- Every row pays it, so the total is ≥ max(β + 3^{1/3}, β + P, 2β).
- The minimum is **2.1144**, at β = 0.672.
- Improved DL initial splittings lower 3^{1/3}, but they can never go below the bottom bound 1.923.

### 3.2 Theorem K: split the rational side into k primes, then descend

**The device.** Accept only splits whose large primes are k primes from one dyadic interval, then descend
each of them (hypD Theorem K). At fixed degree this is the L[1/2] engine.

**At growing degree [model-PROVED].**
- The rational value N^{1/d}·E splits into primes ℓ ≈ N^{1/(dk)}.
  - With k bounded, ℓ = L[2/3].
  - With ℓ = L[1/3], k ~ u, this is just ordinary smoothness.
- Descending ℓ = L[2/3, y] in the same degree-d field uses the 2-dimensional lattice {ℓ | a − bm}.
  - Its vectors have size ≥ √ℓ.
  - The algebraic norm is then ≥ ℓ^{d/2}·N^{1/d} = exp((δy/2)·x·(1+o(1))) = N^{Θ(1)}.
  - That is far above L[2/3] scale, so there is no finite L[1/3] constant.
- Descending inside the main field is therefore viable only for ℓ ≤ L[1/3, O(1)], which is the
  large-prime regime of row 14.
- Codex heurnfs_1 §2 derives the same conclusion from Theorem K's own exponent system: one level forces
  α ≥ 2/5.

### 3.3 Recursion

**The device.** recurse Theorem RA: replace the rigorous shape test by a recursive call of the algorithm.

**Under H [model-PROVED].**
- At L[1/3] every smoothness or cofactor test is already L^{o(1)}.
  - ECM on primes ≤ L[1/3] costs exp(O(x^{1/6}(log x)^{5/6})).
  - Even NFS on an L[2/3]-size norm costs exp(O(x^{2/9}(log x)^{7/9})) = L^{o(1)}.
- In MNFS the test term (1.777) is not binding. Recursion changes the price of *recognising* a
  candidate, never the *probability* that a candidate is a relation.
- Large-prime variations with ≤ k primes ≤ B^λ raise the probability by ρ(u − O(k))/ρ(u) = u^{O(k)}
  = L^{o(1)} [KNOWN; e.g. Lenstra–Manasse, Dodson–Lenstra, Lambert].
- Recursion does change constants at exponent 1/2 in the *rigorous* model (RA: 0.8510 at d = 9), because
  there the shape test is L[1/2, 1/√2] and binds. That is not an L[1/3] effect.

### 3.4 The N·w multiplier

**The device.** ldbound: F(g₀, g₁) = N·w, with w a random prime slightly above the rational values.

**Under H [model-PROVED].**
- *Literal design.* gbe(Nw) leaves the heights of a₁, …, a_d at M and puts w into a₀ ≈ Mw.
  - With w ≈ |G| ≈ N^{1/d}, the algebraic norm gains a factor N^{1/d}, so κ: 2 → 3.
  - GNFS becomes (32/3)^{1/3} = **2.2013**, MNFS **2.1909** (`mnfs(ka=2)`).
- *Re-expanded design.* Expand Nw in base (Nw)^{1/(d+1)} instead (codex): all heights grow by
  w^{1/(d+1)}, and c = 1.9230·(1+ω)^{1/3} for w = N^ω.
  - In ldbound w = N^{1/d+o(1)}, so ω → 0 and the effect is o(1).
- *Many multipliers k_iN sharing m* are just more polynomials with root m mod N. That is MNFS, already
  optimised: V = L[0.125] costs heights a factor V^{1/(d+1)} = (log N)^{O(1)}.

### 3.5 The 2D-random rational side, and multiple rational sides

**2D-random rational side [model-PROVED].** pcest (a″): g = g₁x − g₀, with (g₀, g₁) uniform in [M/2, M]².
The heights are those of base-m, so the ledger and the constant are unchanged (1.9230 / 1.9019). The device
controls pair correlations, which H assumes away.

**Multiple rational sides [PROVED, elementary].**
- Two linear forms a_ix − b_i with the same root mod N satisfy N | a₁b₂ − a₂b₁. If 2H₁H₂ < N they are
  proportional (codex heurnfs_1).
- So at heights N^{1/(d+1)} there is only one rational side up to scaling.
- Other options:
  - A degree-1 side for *another* root of f mod N needs that root, which needs the factorization.
  - Giving f W prescribed roots m₁, …, m_W raises its height to N^{W/(d+1)} and gains only a factor W.
    That is worse for every W ≥ 2 that is O(1), and impossible for W = L[·].
  - Reusing one f with many 2D-random g's needs F(g₀^{(i)}, g₁^{(i)}) ≡ 0 mod N. That has only
    ≈ M²/N ≪ 1 expected solutions at height M.
- MNFS's asymmetry (one cheap shared side, many expensive sides) is therefore forced.
- The symmetric alternative, "any two of V degree-d sides", costs **2.3962** (`mnfs_symmetric`).

### 3.6 Batch smoothness

**The device.** lv6: product and remainder trees replace HECM.

**Under H [model-PROVED].** Testing m values against B costs (m + B)·L^{o(1)}, exactly ECM's per-candidate
L^{o(1)}. The number of smooth values is unchanged, and the test term does not bind (§1). No change.

### 3.7 Norm-height versus search-entropy separation (r20 #5)

See §5. Generic higher-dimensional and toric supports reproduce the 2x/D + D·C ledger, so there is no change. The
special family that r20 needs is not constructed.

## 4. Barriers (Step 3)

### 4.1 Barrier B1: the two-polynomial height LP [model-PROVED under G]

**Setting.** Take polynomials f₁, f₂ of degrees d₁ ≤ d₂ with heights N^{h₁}, N^{h₂} and a common root mod N.

**Constraints.**
- *Resultant [PROVED].* N | Res(f₁, f₂) ≠ 0, and |Res| ≲ H₁^{d₂}H₂^{d₁}. So d₂h₁ + d₁h₂ ≥ 1.
- *Counting [G].* (d₁+1)h₁ + (d₂+1)h₂ ≥ 1.

**The LP.** `height_lp` solves it exactly:
- min(h₁ + h₂) = 2/(d+2) for (1, d);
- min(h₁ + h₂) = 1/d for (d, d);
- min(h₁ + h₂) = 2/13 for (2, 10).

In every case the minimum is 2/D·(1 + O(1/D)), with D = d₁ + d₂.

**Consequence.** The norm product is ≥ N^{2/D}·E^{D}, so log X ≥ 2x/D + D·log E ≥ 2√(2x·log E).
- Codex rightly notes that (R) alone does *not* give a universal lower bound on H₁H₂. The counting
  heuristic G is what does that work.
- Under G every two-polynomial selection is GNFS-equivalent at leading order.
- SNFS shows that G is the operative assumption: dropping the 2x/D to x/D gives 1.526.

### 4.2 Barrier B2: sieving dimension [model-PROVED]

**Setting.** Use sieving polynomials of degree t − 1, so the search entropy is C = t·log A.

**Ledger.** log X ≈ 2(t−1)x/D + D·C/t. Minimising over D gives 2√(2(t−1)xC/t), which is increasing in t.

**Constants.** c(t) = (128(t−1)/(9t))^{1/3}: 1.923, 2.117, 2.201, … → 2.423. The best choice is t = 2.

**Toric supports.** These have the same structure: every extra coordinate pays the full coefficient height
on both sides. This is r20 #5's "generic version is killed".

### 4.3 Barrier B3: tower rings do not help factoring [model-PROVED under G]

**Setting.**
- Base ring Z[ι] with h(ι) of degree η; the ideal 𝔫 = (N, ι − r) has norm N.
- f is defined over Z[ι] with degree d; the rational side is g = x − m(ι).
- In integer-height units H′ = η·h the constraints are:
  - counting: (d+1)H_f′ + H_m′ ≥ 1;
  - resultant: H_f′ + dH_m′ ≥ 1, because Nm(f(m)) ≥ N when f(m) ∈ 𝔫 \ {0}.
- The norm product is H_f′ + H_m′ + (d+1)C/2.

**Result.** η cancels. The ledger is the 1-D ledger with a monic linear g, at value 9/29 for d = 5 for every
η (`tower_ledger`).

**Why exTNFS gains in DL.** For F_{p^n} the constraint is modulo p, with det p^{(…)} and not p^n, while the
resulting norm still has to be ≥ p^n·(…). Z/N has no proper "subfield modulus", so this saving is not
available.

### 4.4 Barrier B4: what binds in MNFS; the two hypothetical levers

**What binds [model-PROVED].**
- At the MNFS optimum, sieve = LA = 1.9019 and the relation count is tight.
- The test term is 1.777 and has slack.
- An optimistic block-angular LA, max(2β, v + 2β_a) in place of 2·max(β, β_a + v), also gives 1.9019.
  So the per-field column structure cannot be exploited.
- Any device that does not reduce E², the per-relation smoothness penalty, or ω leaves c_M unchanged.

**The two levers [HYPOTHETICAL].**
- **(a) Output-sensitive enumeration of rational-smooth (a, b).**
  - If those pairs could be listed in time ≈ their number, E²·L[−r], the MNFS program falls to **1.6386**.
    This is the factory's per-N constant, which is the factory's own mechanism without the amortization.
  - *Why it is hard [model-PROVED].* Special-q, or "s = s₁s₂" constructions, reduce the examined points
    only by known factors Q ≤ E² = L[1/3, O(1)]. Those change u by O(1), a gain of ρ(O(1)) = Θ(1). So
    every lattice or product trick of this kind stays at E²·Θ(1).
  - Smooth values of a − bm form a set of density L[−0.25] in the box. Enumerating them below box cost is
    an open algorithmic problem. It is also exactly the form a positive r20 #5 result would have to take.
- **(b) Sparse LA exponent.**
  - With ω ≤ 1.5, MNFS drops to 1.8597, where the sieve alone binds.
  - GNFS does not move, since e = β is already forced by the sieve.
  - No sub-quadratic dependency finder for sparse GF(2) matrices is known. Wiedemann/Lanczos costs
    C·nnz = C^{2+o(1)}.

### 4.5 Barrier B5: why the exponent is 1/3 [model-PROVED for the ledger]

**Argument.**
- Barriers B1–B3 give log X ≥ 2√(κ′·x·C) for search entropy C.
- With C ~ x^α, X has exponent (1+α)/2.
- Smoothing at L[α] costs exponent (1+α)/2 − α = (1−α)/2.
- That must be ≤ α, so **α ≥ 1/3** (`alpha_floor_generic`).

**The 1/3 is the balance.**
- The degree must grow like x^{1/3} to shrink N^{1/d}.
- The search entropy must then pay d·log E.

**Other floors.**
- One descent level gives the worse floor 2/5: 1 − η ≤ α and (1+η)/2 ≤ 2α (codex; lemmaM (b)).
- Multi-level descent with O(log log N) contracting levels reaches 1/3 in the size model (recurse §3.4),
  but its bottom level is bound by B1.

**Scope.** This is a barrier for algorithms whose relations are smooth values of resultants of a small search
family against fixed polynomials with a root mod N. It is not a lower bound for factoring.

## 5. Structurally new attempts below 1/3 (Step 3), with their precise barriers

**(i) Toric / entropy–height separation (r20 #5).**
- Below 1/3 the ledger would have to be the *aggregate* log X ≪ x/d + (d/r)·C with r ~ x^{1/4} free.
- That gives d ~ x^{1/2}, log X ~ x^{1/2} and α = 1/4 (`alpha_floor_aggregate`).
- *Barrier.* Generic multivariate supports pay mixed-resultant multidegrees. Each search coordinate costs
  coefficient height on both sides, which returns the ledger of B2.
- *Needed.* A family where many coordinates share one height cost, together with:
  - a bounded rational side and bounded denominators;
  - enough relation rank.
- No such family is known. There is no constant to compute: the formal objective max(β + h/(4β), 2β)
  = √h needs a height h that nobody can supply.

**(ii) Small-characteristic elimination lifted through Witt vectors (r20 #6).**
- *Naive lift [PROVED impossible].* A unital map Z/p^k → Z/N needs N | p^k.
- *Corrected lift.* A lift with corrections needs integer carry terms whose aggregate norm is L[1/2, h].
  Generic carries are N^{Θ(1)}, which forces α ≥ 1/2.
- *What the quasi-polynomial DL relies on.* Kleinjung–Wesolowski (KNOWN) use the Frobenius identity
  x^q − x = ∏(x − a) and a known field of definition. Z/N has neither:
  - its characteristics are unknown;
  - powering creates carries.
- *Prime fields.* F_p with p prime is the DL analogue of Z/N. Its best constant is also c_M (MNFS-DL), and
  nobody has improved it since 1993. This is consistent with the barrier above.

**(iii) Lattice / Coppersmith-style constraints.**
- Polynomial selection is already a lattice problem: degree-d polynomials with root m mod N form a lattice of
  determinant N in dimension d + 1.
- The Gaussian heuristic makes the base-m height N^{1/(d+1)} optimal up to O(1) factors. The skewness
  and α-optimisation of practice are L^{o(1)} gains.
- Coppersmith's small-root method needs extra structure, such as known bits of p. For general N it supplies
  no polynomial below the counting bound G.

**(iv) Recursion or descent below 1/3.** Every descent tree ends in bottom relations bound by B1. Recursion only
reprices recognition (§3.3). Neither touches the ledger.

## 6. GPT-6.1-sol blind round (codex/heurnfs_1.out)

**Agreements.** Run before our analysis, it independently:
- reproduced GNFS, MNFS (with identical optimum parameters), the factory (1.6386/2.0069) and the
  Kim–Barbulescu constants;
- reached the verdict "no improvement on c_M, no exponent below 1/3";
- gave the bottom bound for descent;
- gave α ≥ 2/5 for one-level descent and Theorem K;
- gave the (1+ω)^{1/3} multiplier scaling;
- proved the multiple-rational-sides determinant fact.

**Disagreement.** Codex's multiplier row assumes re-expansion of Nw. We add ldbound's literal design
(w in a₀), which costs κ = 3 (2.2013). Both are correct for their designs.

**Our additions.**
- the exact height LP B1;
- tower cancellation B3;
- the t-dimension closed form;
- the symmetric-MNFS and block-LA programs;
- the two hypothetical levers with their constants (1.6386 and 1.8597);
- the 2.1144 random-square-row bound.

## 7. Reviews

**Pending at hand-off.** An Opus model critic and an Opus literature critic were launched but had not reported when this note was committed. No improvement is claimed, so no claim depends on them; their fixes are still to be applied.

## 8. Reproduction

```bash
.venv/bin/python -m src.heurnfs.models          # Table 1 (≈ 2 min)
.venv/bin/python -m pytest -q tests/test_heurnfs.py
```
