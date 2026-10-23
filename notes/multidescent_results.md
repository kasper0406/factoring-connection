# multidescent: conditional goodness at level ≥ 2 and two-level descent — results

Agent `multidescent`, phase 25, 2026-10-23.

Files:
- This note; pre-registration `notes/multidescent_prereg.md` (commit `80d25ee`, code frozen there).
- Code: `src/multidescent/toy.py` (toy fields, special-q boxes, children, constrained resampling; diagnostic),
  `experiments.py` (E1–E3), `analyze.py`, `model.py` (exact size-model and union-bound bookkeeping).
  Tests: `tests/test_multidescent.py`.
- Data: `experiments/raw/multidescent_{e1,e2,e3d3,e3d5}.jsonl`, summaries `experiments/raw/multidescent_*_summary.json`.
- Codex: `codex/multidescent_1.{md,out}` (blind round, launched before any analysis here), `codex/multidescent_2.*`
  (adversarial round, §9).

Labels: **proved**, **proved conditional on X**, **model** (exponent arithmetic under stated premises), **heuristic**,
**experimentally supported** (toy), **open**, **UNVERIFIED**. Everything about the project's algorithm rests on
internal, unrefereed results (Theorem U and its chain; nested; escapeR2).

---

## 0. Verdict

**Plain language.** We did not obtain a rigorous factoring algorithm below L[2/5]. A two-level descent would give
L[3/8], but it still needs analytic estimates that nobody has proved. The new finding is about *which* estimates.
Earlier notes (nested Lemma OD, escapeR2 Hyp LS_C) framed the two-level proof as a joint-smoothness estimate for
one polynomial at ≈ d^{3/2} points at once, more points than it has coefficients. We show that this joint estimate is
not needed: a union bound over the children of each parent needs only conditional statements about **three
geometric points at a time** (the parent and a pair of points in one child's box). The price is precision. Each
child must fail with probability at most ≈ 1/(number of children), so the 3-point correlations must be known to
relative error ≈ D^{−3/8} (D = log N/log log N), a polylogarithmic but demanding precision. It sits at the scale of
the secondary terms in smooth-number asymptotics, so the estimate must be a genuine covariance bound. Its rational
part is a ternary linear equation in smooth numbers: Harper's a + b = c, but with large coefficients, weights and a
conditional normalisation. The algebraic part is of Bombieri–Vinogradov type for a structured, sparse family of
moduli. The level-1 moments themselves also need a mixed-smoothness re-proof (Hyp U′), which is not routine. Fresh
randomness (new fields, multipliers, resampling) does not escape the obstruction. The union bound replaces joint
control by precision; it does not shrink the certificates (Prop R3 and Lemma OD stand as statements about them).

Precisely:
1. **The children's law (Q1).**
   - *Exact form* [proved; codex blind (1)]: the incidence law of children is a ratio of selection-weighted
     divisibility counts (§2.1). E1: algebraic labels follow the "weighted uniform" law ∝ 1/q per available label to
     within ±11% in every sub-bucket; against uniform-on-labels the ratio stays in [0.80, 1.28] [toy].
   - *Rational children transfer* [Prop RC, proved conditional on Hyp U_Q (Theorem U's moment decomposition at the
     child scale) and Lemma CL (averaged upper bound for the children's law, outline repaired after review,
     UNVERIFIED)]: a pool of (log N)^{O(1)} fresh fields, targetwise Markov upgrades of Theorem U's moments, and a
     removal sieve that uses only an UPPER bound on the children's law give "every rational child is pool-good
     outside density 2ε′ per window" at a factor-2 loss in the level-1 count. No joint statement about children is
     needed. This answers codex's blind objection (codex adversarial agrees, conditionally). By R2/R3 this alone
     cannot go below 2/5.
   - *Algebraic children* [proved, structural]: goodness is a property of (label, parent's field). Conditionally on
     the parent's values the field law lives on a fibration: the line g + Z·P₀ in the rational side, and over each g
     the constrained-resampling fiber of rank d − 1 in c, with weights given by the coefficient box (§2.4).
2. **Exponents and what must be proved (Q2).**
   - *Model*: α_j = (j+1)/(3j+2) (2/5, 3/8, 4/11, 5/14, …), d = D^{j/(3j+2)}, D^{α_j} children per node at every
     level (escapeR2 §7.2; codex blind; `model.py`). A polylog per-child loss is affordable as *search work*; as a
     loss of *root density* it would be raised to the power 2k and is not affordable (codex adversarial #7; §3.1).
   - **Prop UB (union-bound reduction) [proved].** If the two-level count's parents satisfy Theorem-U-type level-1
     moments and the children's good-point counts have small conditional relative variance *on average over
     parent–child incidences* (Hyp PC₃, weighted form), then the two-level descendable density is at least a quarter
     of the level-1 Paley–Zygmund ratio. No statement about two children at once is used (§3.2).
   - **Theorem 2L [proved conditional on PC₃, on Hyp U′ (mixed-smoothness re-proof of FM-U/LD*-U/CS with B-smooth
     rational and B₂-smooth algebraic sides, a hybrid restricted BV; OPEN, plausible) and on the internal chain].**
     A Las Vegas algorithm factors every N in expected L_N[3/8, c], c finite (§3.3).
   - What PC₃ needs (§3.4): relative precision O(μ₁/W) (≍ D^{−3/8} under the typical-factorisation count; ≥ D^{−1/2}
     from the deterministic child cap), as a *covariance* statement, for (i) a rational block — three linear forms in
     the two coordinates of g, a ternary equation with coefficients ≤ x^{o(1)} (Harper's Corollary 1, coefficient 1,
     has relative error O(log(u+1)/log y) of exactly this scale, but from main-term bookkeeping, not minor arcs), and (ii) an algebraic block — three evaluations in (w, c), whose
     reduction to BV needs density-relative estimates for an endogenous sparse modulus family (not established).
     The random model predicts covariance O(D^{−2α}) ≪ η [heuristic]. **We prove no α < 2/5, and found no weaker
     provable statement that gives one.**
   - j ≥ 3 levels: the precision needed at depth i is D^{−iα}; from j = 4 it is below the heuristic local-factor and
     archimedean covariance D^{−2α}, and conditioning on local data cannot repair this without pinning the field
     (§3.5). Deeper descents, and α → 1/3, remain open by this route.
3. **Fresh randomness (Q3)** [§4]: a fresh field per level refreshes rational children only (Prop R; R2/R3 then cap at
   2/5); per-child fields are blocked by Lemma FB/Cor EX; a fiber field per parent gives fresh randomness only to the
   ≤ d − 1 parents it retains (other rows in it are ordinary one-field rows), and private fields fall foul of the
   singleton-column Lemma CC; multipliers w add one coordinate (Prop NB); constrained resampling spends the same
   budget. The rational side is the tightest budget: rank ≤ 2, and 1 after conditioning on the parent. **The
   two-level case escapes the high-point joint estimate only through the union bound (precision), not through
   randomness.**
4. **Barrier (Q4)** [§5]: for the pooled design class the necessary-and-sufficient condition is escapeR2's T_q
   (codex: equivalent to pooled coverage, up to constants). Two sufficient analytic routes are known for 3/8: LS_{3/2}
   (complete-certificate moments, L^{O(1)} ratio, ≈ d^{3/2} dependent forms) and PC₃ (3-point conditional covariance,
   precision ≍ D^{−3/8}). Neither implies the other. Lemma MP makes PC₃'s precision necessary for any moment-based
   union-bound proof. We know no derivation of PC₃'s rational block from level-of-distribution hypotheses (EH
   included); black-box AP statistics do not determine ternary counts for general sets. Its prototype is Harper's
   unconditional a + b = c asymptotic; the missing piece is a conditioned, weighted, growing-coefficient covariance
   version beating the secondary-term scale [open]. Codex's estimate that PC₃ is provable by known methods in a
   focused effort: ≈ 15%; ours: 5–10%.
5. **Experiments** (pre-registered; toy, 48-bit N; 10 predictions, 8 held): E1 children's law ≈ 1/q-weighted
   uniform on labels (P1a, P1b held; P1c failed at the top rational window, a finite-size edge, explained); E2 rational
   children are exactly as good as log-matched uniform primes in fresh pools, also in the lower tail (P2a, P2b held:
   +0.003 [−0.002, 0.008]); E3 conditional child goodness factorises up to field quality, including beyond the fiber
   rank (P3c failed against a poor null; post hoc within height classes the ratio is 0.992), and child-box counts are
   nearly Poisson under the fiber and line laws (P3a, b, d, e held).

---

## 1. Setting, and what is already known

Notation as lv5unif §1 and nested §8: t = ln N, D = t/ln t, Λ_α = D^α ln t, L[α, c] = e^{(c+o(1))Λ_α}. A *field* is
A = (g, w, c) (design (a‴)): G = g₁x − g₀y, F = F̂_{g,Nw} + G·R_c, F(g₀, g₁) = Nw, degree d. A *rational target* ℓ
has lattice Λ_ℓ(g) = {P : ℓ | G(P)}; an *algebraic target* 𝔮 = (q, r) has lattice {P : a ≡ rb (mod q)}. The
children of a box point P are the primes > B of its rational cofactor and the degree-one prime ideals > B of F(P).
An algebraic child of P = (a₀, b₀) has label (q, a₀b₀^{−1} mod q); its lattice is Λ_𝔮 = {P′ : det(P₀, P′) ≡ 0 mod q},
which depends only on (P₀, q).

Prior results this note builds on (not redone):
- **Theorem U** (lv5unif): rigorous L[2/5] via one-level descent; the floor comes from the split acceptance
  (log N)^{−2k} (k ≲ D^α) and the norm exponent θ ≳ 1/√k.
- **nested**: Theorem NR (ERH + AD + LV6′ ⇒ L[1/3]); Prop SI (density alone cannot help in a fixed field); Lemma FB;
  Lemma RW (row algebra for multi-level trees); Lemma OD (a degree-d form has d + 1 coefficients); Prop R2 (no
  algebraic children ⇒ 2/5 in the model).
- **escapeR2**: Cor EX (no export of algebraic children between unrelated fields); Prop R3 (O(d)-point algebraic
  components ⇒ 2/5, model); the chain sub-model (3/8 at component size d^{3/2}); **Hyp LS_C** (moments of complete
  certificates) and **Theorem NR-q** (pooled fields; hypothesis T_q, field-averaged and target-pointwise, success
  L^{−O(1)} on a polylog-dense set of split primes).
- **almostall13**: Prop NB (rank of (F(P_i), G(P_i)) ≤ min(s, d + 1) + min(s, 2)); Lemma MP (moments accurate to
  relative error η cannot certify Pr[dead] < η).
- **navigate §5.4**: the local-backtracking criterion Σ m_i ln(1/γ_i) = O(Λ) and the request for a conditional
  goodness bound "uniform over a set of q of full measure under the children's law".

*Reconciliation with navigate §5.4.* escapeR2's NR-q already shows that a *field-averaged, target-pointwise*
statement suffices at the root, so "uniform over full measure" is not needed for the split primes. Inside a field,
navigate's local-backtracking criterion is one sufficient strategy; §3 shows that its "full measure" form (per-child
failure ≤ 1/m) is exactly what makes a 3-point argument possible.

## 2. The children's law (Q1)

### 2.1 Exact form [proved; codex multidescent_1 (1)]

Fix the field and the history h. Let b_h(P) ≥ 0 encode parent selection (acceptance, first-success rules) and
n_τ(P) the multiplicity of target τ among P's children. The incidence law of children in a window J is
ν_h(τ) = Σ_P b_h(P)n_τ(P) / Σ_P b_h(P)Σ_{τ′∈J}n_{τ′}(P). Choosing one random child per parent adds the factor
1/m(P). For the *removal sieve* below only the counting measure (all accepted parents, all their children) matters,
and only through an upper bound.

**Heuristic shape.** For a random value n of size X conditioned to be y-smooth, Pr[q | n] ≈ Ψ(X/q, y)/Ψ(X, y) ≈ q^{−β}
(β = β(X, y), Hildebrand–Tenenbaum), and an algebraic label (q, r) needs P ≡ r (mod q) projectively (probability
1/(q + 1)). So the label weight is ≈ q^{−1}·q^{1−β}, which varies by ≤ 2^{1−β} + o(1) inside a dyadic window: the
children's law is "uniform on available labels, tilted by 1/q". This breaks down only when q is comparable to the
value (X/q small), where Ψ(X/q, y) collapses — the rational children of E1's top window (§6.1).

**Codex's warning (4) [proved there].** In a *fixed* field the children of a search with W candidates lie in a set of
≤ 2W log X/log B labels, so no polylog domination by the uniform law holds pointwise in the field. All domination
statements below are averaged over the field (and, for rational children, the pool is independent of the field).

### 2.2 Domination of the children's law by the uniform law [outline, UNVERIFIED; experimentally supported]

*Design RC.* Theorem U's design with the level-1 rational cofactor allowed to be B₂-smooth, its prime factors in
(B, B₂] being the rational children. (In Theorem U's own design the rational side is B-smooth and the rational
statement below is vacuous; proof critic M1.)

**Lemma CL (averaged upper bound).** In design RC, for every prime q ∈ (B, B₂],
E_{ℓ,A}[#{P : Acc_ℓ(P), q | G(P)/ℓ}] ≤ (log N)^{C}·q^{−β₂}·E_{ℓ,A}[#{P : Acc_ℓ(P)}], β₂ := β(MV/Y, B₂).

*Proof outline (repaired after the proof critic's M1 and codex adversarial #6).* The first draft bounded the
algebraic probability at good moduli by `def:good` and discarded bad moduli by Lemma 8.5′. That is invalid for an upper
bound: Lemma 8.5′ bounds the *number* of bad rational candidates, with relative error ε_N = exp(−cu/log²u) + y^{−1/2},
which is far larger than ρ_F, so a bad modulus with algebraic probability up to 1 is not negligible. Repair: use the
w-averaged bound valid at every modulus, as in LD*-U Steps 2–4. Lemma W-U gives Pr_w[w ≡ r (mod t)] ≤ C_W/φ(t), so
Pr_{w,c}[F(P) acceptable] ≤ (log N)^{C}ρ_F for every P (DS-U/L2-U-type counting of smooth values in the layer). On the
rational side, FIB (|𝒢_P(n)| ≤ 2M/V) and the number of B₂-smooth multiples of q in I_P, ≤ Ψ(z₂/q, B₂) ≤ Cq^{−β₂}Ψ(z₂, B₂)
(Hildebrand–Tenenbaum Theorem 3, 1 ≤ q ≤ B₂), compared with the good-mass lower bound cλρ of [pp] Thm FM. Constants are
not tracked; the divisor-weighted bad-mass issue codex raises is avoided because no bad-modulus subtraction occurs
in an upper bound of this form. The algebraic analogue (q | F(P)) needs the residue mod q from c combined with W-U
mod t by CRT; not written. E1 measures the toy domination constant: ≤ 1.28 for algebraic labels.

### 2.3 Prop RC: rational children transfer, with explicit loss [proved conditional on Hyp U_Q]

Rational children are integers, so each can be descended in a fresh field (nested Prop R; multi-field rows as in
lv5unif §11/escapeR2 §8).

**Hyp U_Q (Theorem U's structure at child scale Q).** For targets q in a dyadic window (Q, 2Q] ⊂ (B, B₂] and the
one-level count N′_q(A′) in a fresh field A′: (a) E_{A′}N′_q ≥ 𝔤 − 𝔟_q with 𝔤 > 0 independent of q and
(1/#window)Σ_q 𝔟_q ≤ ε_N𝔤, ε_N = exp(−cu/log²u) + y^{−1/2}; (b) (1/#window)Σ_q E_{A′}[N′_q(N′_q − 1)] ≤ K𝔤². This is
exactly the form of [pp] Thm FM (good mass pointwise in the target, bad mass bounded on average by Lemma 8.5′) and
Thm LD*. (Re-checking (P1)–(P10) at scale Q is routine but UNVERIFIED; for Q < B^{K} the factor base of the pool
fields must be lowered to Q^{1/K}.)

**Prop RC.** Work in design RC. Assume U_Q for every window with K ≥ 1 and K + 𝔤^{−2} ≤ (log N)^{O(1)}, and Lemma CL.
Let ε′ ∈ (4ε_N, 1). Then:
1. *(targetwise upgrade)* every q outside a set of density ≤ ε′ in its window has
   p_q := Pr_{A′}[N′_q ≥ 1] ≥ ε′/(32(K + 𝔤^{−2}));
2. *(pool)* with J_f := ⌈32(K + 𝔤^{−2})ε′^{−1}·ln(8t/ε′)⌉ independent pool fields, with probability ≥ 1/2 every
   window has at most a 2ε′ fraction of pool-bad targets;
3. *(removal sieve)* E[#{accepted parents all of whose rational children are pool-good}] ≥ E N₁·(1 − 2ε′(log N)^{C_ν}),
   where (log N)^{C_ν} := (log N)^{C}·Σ_{windows Q} Q^{1−β₂}/ln Q is polylogarithmic (Lemma CL's constant summed
   over the ≤ 2t windows; Q^{1−β₂} ≤ B₂^{1−β₂} ≍ u log u by Hildebrand–Tenenbaum). With ε′ := (log N)^{−C_ν}/4 the
   level-1 count loses at most a factor 2.

*Proof.* (1) By Markov, #{q : 𝔟_q > 𝔤/2} ≤ 2ε_N·#window ≤ (ε′/2)#window, and #{q : E[N′_q(N′_q−1)] > 2K𝔤²/ε′} ≤
(ε′/2)#window. For the remaining q, E N′_q ≥ 𝔤/2 and E[N′_q(N′_q − 1)] ≤ 2K𝔤²/ε′, so by Paley–Zygmund
p_q ≥ (E N′)²/(E[N′(N′−1)] + E N′); if E N′ ≤ 1 this is ≥ (𝔤²/4)/(2K𝔤²/ε′ + 1) ≥ ε′/(8(K + 𝔤^{−2})), and if E N′ ≥ 1
it is ≥ (E N′)²/(2K𝔤²/ε′ + E N′) ≥ min(1/2, ε′/(16K)). (2) For such q, Pr[q pool-bad] ≤ (1 − p_q)^{J_f} ≤ ε′/(8t);
the expected pool-bad fraction of good q in each window is ≤ ε′/(8t), and Markov plus a union bound over the ≤ 2t
windows gives probability ≥ 3/4 that every window has ≤ ε′ further bad targets (so ≤ 2ε′ in all). (3) Every accepted parent with a
pool-bad rational child is counted once per such child; by the averaged domination, the expected number of
(parent, child) incidences with the child in a set S of density ≤ 2ε′ of the window (Q, 2Q] is, by Lemma CL,
≤ Σ_{q∈S}(log N)^{C}q^{−β₂}E N₁ ≤ 2ε′(log N)^{C}Q^{1−β₂}(ln Q)^{−1}·2E N₁. Sum over windows. The pool is independent of
(ℓ, A), so the pool-bad set may be fixed before averaging. ∎ (Step 1's constants: if E N′ ≤ 1 the ratio is
≥ (𝔤²/4)/(2K𝔤²/ε′ + 1) ≥ ε′/(8(K + 𝔤^{−2})); if E N′ ≥ 1 it is ≥ min(1/2, ε′/(16K)); the stated ε′/(32(K + 𝔤^{−2})) is
weaker than both.)

**Explicit loss.** Count × 1/2; pool size (log N)^{O(1)}; per-child goodness from ε′(log N)^{−C} to "all but density
2ε′" per window; ε′ = (log N)^{−C_ν}/4. Only an UPPER bound on the children's law is used. *Why codex's objection does not apply:*
Theorem U's γ is a joint statement, but its proof gives the good mass *pointwise in the target* and the bad mass on
average, and the second moment is a non-negative average, so Markov makes the exceptional set as small as any
polylog at a polylog price in p_q, which the pool then removes. Codex's statement (2) (domination gives only
ν(Bad) ≤ K·U(Bad)) is exactly what step 3 uses. Codex's adversarial round agrees that this answers the blind
objection logically, conditional on U_Q at every child scale and on Lemma CL.

**Scope.** In the 3/8 chain the level-1 rational sides are B-smooth anyway (rational value exponent 1 − δ = 2α), and
by R2/R3 rational children cannot by themselves move α. Prop RC answers Q1 for rational children; the exponent is
decided by the algebraic children.

### 2.4 Algebraic children: the conditional law [proved, structural]

Fix ℓ and a level-1 point P₀ = (a₀, b₀), and condition on the parent's exact values G(P₀) = n_r, F(P₀) = n_a.
- *Rational side.* {g : g₁a₀ − g₀b₀ = n_r} is the line g⁰ + Z·P₀. Along it G(P′) = G⁰(P′) + j·det(P₀, P′), and for
  P′ ∈ Λ_𝔮, q | det(P₀, P′). So each child point's rational value runs through an arithmetic progression in j whose
  modulus is a multiple of q (one free direction).
- *Algebraic side.* For fixed (g, w), {c : F(P₀) = n_a} is an affine lattice of rank d − 1 (when nonempty), i.e. the
  constrained-resampling fiber F + G·L₀·H (nested §4; L₀ = b₀x − a₀y): every child label (q, a₀/b₀) survives and stays
  a root of every fiber polynomial (E3: 0 violations).
- *This is a fibration, not a product* (codex adversarial #3; proof critic minor 4): F̂_{g,Nw} depends on g, the
  c-coset over a given (g, w) is nonempty only when G(P₀) | n_a − F̂_{g,Nw}(P₀), and the conditional weight of (g, w) is
  proportional to the number of compatible c in the coefficient box. Every analytic statement below uses the literal
  conditional law.
- *Consequence.* A child's goodness is a function of (label, field); its conditional probability given the parent is
  a statement about this conditional law. Transferring the *label* law (§2.2) is necessary but not sufficient. The
  question is whether one can control the children's goodness *jointly* (all m of them) under this law; §3 shows that
  jointly is not needed, only per child, with precision.

## 3. Exponents and what must be proved (Q2)

### 3.1 The size model [model; `model.py`, tests]

j-level chain (escapeR2 §7.2 with r = j − 1; codex blind derived the same): level-i targets L[ζ_i], ζ₀ = 1 − α
(split), ζ_j = α (factor base); d = D^δ; base-m heights (exponent 1 − δ); smoothness of a value of exponent v over
L[ζ′] at cost L[α] needs v − ζ′ ≤ α. Constraints: 1 − δ ≤ 2α and δ + ζ_i − ζ_{i+1} ≤ α. Summing, jδ + 1 − 2α ≤ jα,
so **α_j = (j+1)/(3j+2)**, δ_j = j/(3j+2), ζ_i = 1 − α − i(3α − 1). Every level has D^{δ + 3α − 1} = D^{α} large algebraic
children per node. Two levels: α = 3/8, d = D^{1/4}, ζ₁ = 1/2, level-1 algebraic value D^{7/8}, level-2 values
D^{3/4}, m ≍ D^{3/8} = d^{3/2} children. The model bounds the exponent of designs of this shape; it proves no supply.

*Per-child losses (corrected after codex adversarial #7).* A factor (log N)^{−Cm} from requiring all m children to
succeed is affordable as **search work per target** (retry until a parent with good children is found:
exp(O(C D^α ln t)) = L[α, O(C)]). It is **not** affordable as a loss of **root density**: the split raises the
descendable density to the power 2k ≍ D^{α}, so (log N)^{−Cm·2k} = exp(−Θ(D^{2α} ln t)) is outside L[α]. navigate's
bookkeeping implicitly assumes the first reading, with independent children (nested Prop BT's model). Prop UB
avoids the issue: under PC₃ the two-level count keeps a constant fraction of the level-1 first moment.

### 3.2 Prop UB: the union-bound reduction [proved]

*Design 2L* (two levels at α = 3/8; the precise parameters are those of lv5unif §1 with d ≍ D^{1/4}, k ≍ D^{3/8},
B = e^{aΛ}, an extra child bound B₂ = e^{b₂D^{1/2}ln t}, and a child-box radius S₂). For ℓ ∈ J and a field A:
- Acc_ℓ(P): P ∈ 𝒫 ∩ Λ_ℓ(g), G(P)/ℓ ∈ 𝒮_B, and F(P) = s·q₁⋯q_m with s B-smooth and q_i distinct primes in (B, B₂] not
  dividing b·D_A (children 𝔮_i = (q_i, ab^{−1} mod q_i)); m(P) := m.
- Box(𝔮): the primitive points, modulo ±1, of the reduced-basis box of Λ_𝔮 of radius S₂ (primitivity matters:
  Good(jP′) = Good(P′) for small j would add perfectly correlated pairs; proof critic minor 1).
- Good(P′) for P′ ∈ Box(𝔮): q ∥ F(P′), G(P′) ∈ 𝒮_B, F(P′)/q ∈ 𝒮_B, plus the gcd checks.
- N₁(ℓ, A) := #{P : Acc_ℓ(P)}, N_𝔮(A) := #{P′ ∈ Box(𝔮) : Good(P′)},
  N₂(ℓ, A) := #{P : Acc_ℓ(P) and N_𝔮(A) ≥ 1 for every child 𝔮 of P}.
- For a triple (ℓ, P, q) let E = E_{ℓ,P,q} := {A : Acc_ℓ(P) holds and q is a child of P} (an event in the field law,
  ℓ fixed). Probabilities are over ℓ uniform in J and A ∼ 𝒜.

**Hypothesis PC₃ (weighted form; suggested by codex adversarial #4).**
  Σ_{(ℓ,P,q)} Pr[ℓ, E_{ℓ,P,q}]·ϕ_{ℓ,P,q} ≤ μ₁/2,  where ϕ := Var(N_𝔮 | E)/E[N_𝔮 | E]² if E[N_𝔮 | E] > 0, and ϕ := 1 otherwise,
and μ₁ := E N₁. Here Var(N_𝔮|E)/E[N_𝔮|E]² = 1/E[N_𝔮|E] + (E[N_𝔮(N_𝔮−1)|E] − E[N_𝔮|E]²)/E[N_𝔮|E]², and
E[N_𝔮(N_𝔮 − 1)|E] = Σ_{P′≠P″} Pr[Good(P′) ∧ Good(P″) | E] concerns the events E ∧ Good(P′) ∧ Good(P″).

*A sufficient pointwise form, PC₃(η, ε):* outside a set 𝔈 of triples with Σ_𝔈 Pr[ℓ, E] ≤ εμ₁, (a) E[N_𝔮|E] ≥ 4/η and
(b) E[N_𝔮(N_𝔮 − 1)|E] ≤ (1 + η/2)E[N_𝔮|E]²; then ϕ ≤ 3η/4 there, and PC₃ follows when (3η/4)W + εμ₁ ≤ μ₁/2. Codex notes
that the pointwise form forces an exceptional incidence fraction ≍ 1/m, so the weighted form is the natural one.

**Proposition UB.** Let μ₁ := E N₁ > 0 and W := E[Σ_{P : Acc} m(P)] (expected number of parent–child incidences;
Σ_{(ℓ,P,q)} Pr[ℓ, E_{ℓ,P,q}] = W). If PC₃ holds, then E N₂ ≥ μ₁/2 and
  Pr_{ℓ,A}[N₂ ≥ 1] ≥ μ₁²/(4·E[N₁²]).

*Proof.* Pointwise, N₂ ≥ N₁ − Σ_{P}Σ_{𝔮 child of P} 1[N_𝔮 = 0] (an accepted parent that is not counted in N₂ has at least
one dead child). Taking expectations, E N₂ ≥ μ₁ − Σ_{(ℓ,P,q)} Pr[ℓ, E ∧ N_𝔮 = 0]. For each triple, Chebyshev in the
conditional law gives Pr[N_𝔮 = 0 | E] ≤ ϕ_{ℓ,P,q} (trivially if the conditional mean is 0). Hence
Σ Pr[ℓ, E ∧ N_𝔮 = 0] ≤ Σ Pr[ℓ, E]ϕ ≤ μ₁/2. Finally 0 ≤ N₂ ≤ N₁, so by Paley–Zygmund
Pr[N₂ ≥ 1] ≥ (E N₂)²/E[N₂²] ≥ (μ₁/2)²/E[N₁²]. ∎

*Well-definedness* (both reviewers): for fixed (ℓ, P, q), E_{ℓ,P,q} is an event in the field law; membership
P ∈ Λ_ℓ(g) is part of it, and on it the child label (q, ab^{−1} mod q), hence Box(𝔮), is determined by (P, q).

**Remarks.**
- No statement about two children at once is used: their failure events are *added*, not intersected. The
  configurations have three **geometric points** (P, P′, P″) and up to six evaluations (G and F at each), plus the
  global acceptance conditions (discriminant, excluded primes). The rational part is three linear forms in the two
  coordinates of g; the algebraic part three evaluations in (w, c), of full rank in c for d ≥ 3 and pairwise
  non-proportional points (Kern-U's kernel argument with Q₃ = L_P L_{P′} L_{P″}).
- The price is precision: the weighted average of ϕ must be ≤ μ₁/(2W). Under the typical-factorisation count
  W/μ₁ ≍ D^{3/8} (so precision ≍ D^{−3/8}); the deterministic cap m(P) ≤ m̄ = O(D^{1/2}) only guarantees that
  precision D^{−1/2} suffices (codex adversarial #5). By Lemma MP (almostall13) any moment-based proof that a child
  fails with probability < η needs moments accurate to relative error η, so this precision is necessary for the
  route. It is polylogarithmic, not L-scale: the union is over the children of one parent, not over all e^{Θ(Q)}
  targets of a field (AD's union bound, nested §6 / almostall13, which needs L[2/3] precision).
- *Precedent.* Removal arguments that sum bad incidences have precedent in rigorous descent: GKZ's proof of their
  Theorem 3 removes the ≤ q^{kd/2+4} values producing trap roots from ≥ q^{kd−4} solutions; KW avoid a bounded-degree
  trap subvariety (codex adversarial #11; literature critic). Those counts are deterministic and hold pointwise for every
  good target — the AD-type property unavailable here — so the analogy is formal. We found no prior use of a conditional 3-point smooth-number reduction of this kind.
- E[N₂²] ≤ E[N₁²] means the second moment needed is the *one-level* one (LD*-type, upper bound with polylog loss).

### 3.3 Theorem 2L [proved conditional on PC₃, U′ and the internal chain]

**Hyp U′ [OPEN; expected but not routine — both reviewers].** Theorem U's FM-U, LD*-U, D-U hold for Acc_ℓ at the
3/8 parameters, with computable box radii S₁, S₂, giving μ₁ ≥ (log N)^{−C}, E[N₁²] ≤ (log N)^{C}(μ₁² + μ₁), search
sizes W₁, W₂ = L[3/8, O(1)] and E[N_𝔮|E] large enough for PC₃. What has to be re-proved:
- *Mixed smoothness.* The rational cofactor is B-smooth, the algebraic value B₂-smooth. Lemma 8.5′(R′) uses one y for
  the values and for the smooth parts of the restricted moduli. With y = B₂ it controls a much larger modulus family
  whose error e^{−cD^{3/8}/log²D} does not beat ϱ(z, B₂)/ϱ(z, B) = e^{Θ(D^{3/8} log D)}; with y = B the algebraic
  u ≍ D^{1/2} violates (R′3) (u log(4u log B)/log B ≍ D^{1/8} → ∞) (codex adversarial #1). A hybrid restricted BV —
  B-smooth modulus cofactors, B₂-smooth values, normalised by the density of the actual modulus family — is needed;
  GM-U, L2-U, DS-U must be redone with two smoothness parameters.
- *Box sizes.* Lemma CS evaluates both Rankin factors at B; at the root algebraic size D^{7/8} that gives an
  L[1/2]-scale box. A mixed (B, B₂) scan replacing CS, and an explicit child radius S₂, are needed (codex #2).
- *Squarefree large part.* The lower bound must subtract values with p² | F(P), p ∈ (B, B₂]: either smooth counts in APs
  mod t·p² (not given by `def:good`) or a union bound Pr_c[p² | F(P)] ≲ 1/p², usable only when 1/B ≪ ρ_F(log N)^{−C}
  (a constraint between the constants; proof critic M3).
- *Value factorisations* up to B₂ = L[1/2] must be recursive ([rec] Lemmas W/OS; the bounded-depth search recurrence
  of escapeR2 §7.3), and Lemma IR must be re-checked at γ₂.

**Theorem 2L.** Assume U′ and PC₃. Then FACTOR_2L —
FACTOR_U with the descendable set 𝒟_A := {ℓ : N₂(ℓ, A) ≥ 1}, DESCEND_ONCE replaced by the bounded-depth search below,
and POSTPROCESS on two-level rows — is Las Vegas and factors every N in expected time L_N[3/8, c] for a finite c.

*Search.* Scan the level-1 box of ℓ; for each P with Acc_ℓ(P), scan the box of each child; return the first P all of
whose children have a good point, with those points. Given U′'s search sizes, the cost is ≤ W₁(1 + m̄W₂) box points
(m̄ = O(D^{1/2})), each with recursive factorisations of values ≤ N^{o(1)} (escapeR2 §7.3's recurrence): L[3/8, O(1)].

*Proof sketch.* (i) Prop UB and U′ give γ₂ := Pr[N₂ ≥ 1] ≥ (log N)^{−C}. (ii) Theorem K-U and Prop select-U hold for
an arbitrary field-dependent target set (Lemma S″), so with δ_A := #𝒟_A/#J they apply verbatim with γ₂ in place of
γ; split acceptance (γ₂/(16 log N))^{2k}/(32k²) = L[3/8, −O(κ)]. (iii) Rows: nested Lemma RW (each algebraic child is
cancelled by its level-2 element z_{P′}; rational sides are B-smooth; |exponents| ≤ (C log N)^3). (iv) LV6 at d ≍ D^{1/4}:
lv5unif Prop LV6-U (factorial-bit characters, cost e^{O(d log d)} = e^{o(Λ_{3/8})}); character primes exceed B₂.
(v) Fruitfulness: Theorem DD with nested Lemma T_multi (the search is deterministic given canonical factorisations).
(vi) Recursion: halves ≤ √(2N) and values ≤ N^{θ+o(1)}, θ → 0; lv5unif Lemma F_α at α = 3/8 suffices for finiteness,
and escapeR2 NR-q step 8's three-term recursion gives a better constant. ∎

**Dependencies.** Internal: Theorem U's chain (lv5unif §10 list, including [pp] Thm FM/LD*, Lemma 8.5′, Kern-U,
L2-U, Theorem DD, Lemma S″, Lemmas OS/W/F/T_rec, Lemma IR), Prop LV6-U, nested Lemmas RW, TS, T_multi, escapeR2 §7.3
(bounded-depth search recurrence); this note's Prop UB. Hypotheses: PC₃ (open), U′ (open; hybrid BV, mixed CS, S₂,
squarefree part). Published: as Theorem U.

### 3.4 What PC₃ asks [analysis; heuristic where marked]

Write the 3-point event E ∧ Good(P′) ∧ Good(P″) as a rational block and an algebraic block.
- **Rational block** (functions of g only): G(P)/ℓ ∈ 𝒮_B, G(P′) ∈ 𝒮_B, G(P″) ∈ 𝒮_B, with the Plücker identity
  det(P′, P″)·G(P) + det(P″, P)·G(P′) + det(P, P′)·G(P″) = 0, i.e. ℓ·det(P′, P″)·x₀ + det(P″, P)·x₁ + det(P, P′)·x₂ = 0
  for the smooth variables x₀ = G(P)/ℓ, x₁, x₂ (codex #11). All three minors are divisible by q; |det(P′, P″)| ≤ 2V₂² and
  |det(P, P′)| ≤ 2V₁V₂ (upper bounds; not uniform two-sided sizes). At 3/8 the values have size e^{Θ(D^{3/4}ln t)} and
  the coefficients, after removing common factors, at most e^{O(D^{5/8}ln t)}: **a ternary linear equation in smooth
  numbers with coefficients ≤ x^{O(D^{−1/8})} = x^{o(1)}**, counted over g in a box (together with the image lattice
  of g ↦ (x₀, x₁) and the transformed box), with conditional weights, to covariance precision O(η). Conditioning on
  G(P) leaves one direction: two linear polynomials in one variable (smooth pairs on a line).
  - Prototype [literature, checked by codex adversarial]: Harper (Compositio 2016, Corollary 1) gives
    #{a + b = c, all y-smooth ≤ x} = (Ψ(x, y)³/2x)(1 + O(log(u+1)/log y)) for (log x)^K ≤ y ≤ x. At this design's scale
    O(log(u+1)/log y) = O(D^{−3/8}), exactly the precision needed. But (literature critic) Harper's minor arcs
    contribute only ≪ (Ψ³/x)(log x)^{−4}; the log(u+1)/log y comes from simplifying the major-arc main term
    (replacing V(x, y; q, δ) by (Ψ/x)V(x, x; q, δ)), i.e. it is main-term bookkeeping that a covariance statement
    cancels. So for coefficient 1 the precision is not the obstruction; the genuinely new difficulties are coefficients
    up to x^{o(1)} (far beyond Harper's major-arc cutoff log^{20} x, so dilations scramble the arcs), the conditional
    weights, and the family average. Under GRH, Lagarias–Soundararajan (Proc. LMS 2012, Thm 2.1) need
    (log x)^{8+δ} ≤ y ≤ e^{(log x)^{1/2−δ}} and fixed smooth weights (their introduction has a one-sentence remark that
    fixed coefficients aX + bY + cZ = 0 should work); Drappeau (MPCPS 154 (2013), Thm 1.3, GRH) removes the upper
    limit, giving relative error O(1/u) on a range that contains y = L[3/8] — again exactly our scale, again of
    saddle-point (main-term) nature.
  - **Hyp RB (rational block) [OPEN].** For the 2L configuration family, with the conditional weights, the 3-form count
    divided by the product of the corresponding 1- and 2-form counts is 1 + O(D^{−3/8−ε}) on (weighted) average.
- **Algebraic block** (given g): F(P) (q | F(P), B₂-smooth with simple children), F(P′)/q, F(P″)/q ∈ 𝒮_B. The map
  c ↦ (R_c(P), R_c(P′), R_c(P″)) has kernel Q₃·Z[x, y]_{d−4} (Gauss's lemma, as Kern-U) and an image lattice whose
  structure is fixed by q and the minors. Peeling one coordinate at a time turns the count into sums of smooth numbers
  in APs with moduli G(·)·(minor factors) and residues determined by the earlier coordinates, at level < 1/2 if H is
  raised by a constant factor in the exponent. **This reduction is not established** (codex #4): full rank gives an
  image lattice, not the needed weighted modulus distribution, reduced residues or multiplicity control. Ordinary BV
  sums over all moduli, whereas this family is endogenous and exponentially sparse, so a density-relative estimate is
  needed; Harper's maximal-residue BV (arXiv 1208.5992, Theorem 1; level √Ψ(x, y)) carries a large-sieve term
  √Ψ·Q·log^{7/2}x, not only the exponential error; Drappeau's Lemme 5 (Compositio 2015) shows that term comes only from
  conductors > x^η, the natural split for a density-relative version (still an absolute bound); Pascadi's level
  5/8 − ε is for a fixed residue, Drappeau's Théorème 1 reaches 3/5 − ε for residues a₁ā₂ with |a_i| ≤ x^δ. [open; plausibly approachable by Theorem-U-type
  arguments with lower and upper bounds]
- **Main terms** [heuristic]: under the random model with local factors, the pair correlation exceeds the product of
  marginals through (i) shared small-prime congruences of pairs (P′ ≡ λP″ mod p, frequency 1/p), (ii) the fields' local
  data (Murphy-α type) and (iii) archimedean size fluctuations (|g|, |c|, magnitudes of the values). Each changes a
  smoothness probability by a factor exp(O(ξ(u)/log y)) = 1 + O(D^{−α}), so if these perturbations have bounded variance
  *uniformly under the parent conditioning*, the covariance is O(D^{−2α}) = O(D^{−3/4}) ≪ η. That uniform bounded
  variance is itself an assertion, not a consequence of the model (codex #5).

**A precision caveat [analysis].** At 3/8 the needed relative precision η ≍ D^{−3/8} is of the same order as
1/u ≍ log u/log y, the typical size of *secondary terms* in smooth-number asymptotics (saddle-point corrections,
Hildebrand–Tenenbaum). So PC₃ cannot be obtained by proving "count = main term × (1 + O(log u/log y))" for the
3-point count and the 2-point counts separately: the constants would have to be compared. It has to be proved as a
*covariance* statement, in which the common main-term shape cancels. In circle-method language: the principal major
arc reproduces the product of the marginals, and what must be shown is that non-principal major arcs and minor arcs
contribute o(η) relatively. (For the algebraic block the analogous statement is that the BV-type error terms, of
size exp(−cu/log²u) + y^{−c}, are what remain once both sides are expanded around the same main term.)

**What Prop UB does and does not avoid** (both reviewers). A union-bound certificate still has 1 + m ≍ d^{3/2} points;
R3 is a model statement about the size of the algebraic component, independent of how its existence is proved, and
it is not affected (R3's O(d) premise is simply not met). Lemma OD's linear relations among those values remain true.
What UB avoids is the *proof obligation* of a joint-distribution estimate at d^{3/2} points (LS_C, or a first moment
of complete certificates): the joint event is never computed, only per-child conditional failure probabilities. In
exchange it needs PC₃-level precision (Lemma MP).

### 3.5 j levels [model + heuristic; corrected after review]

The argument applies level by level: a child at depth i is good if its box has an accepted point all of whose
children are good. With N^{full} ≤ N^{acc}, Var(N^{full})/E[N^{full}]² ≤ (1 + η′)/(1 − ε)² − 1 ≈ η′ + 2ε with ε ≤ m·η_{next},
so depth i needs precision η_i ≲ 1/∏_{i′≤i} m_{i′} ≍ D^{−iα}, and its configurations are the ancestor chain plus a pair
(i + 2 geometric points; the rational block has i + 2 forms in 2 variables). Against the heuristic covariance D^{−2α}
(local data and archimedean fluctuations) the margin (`model.union_bound_margin`) is α, 0, −α, … at depths 1, 2,
3, …: positive for j = 2 (3/8), borderline (constants decide) for j = 3 (4/11), negative from j ≥ 4.

The draft proposed conditioning on the fields' local data (f mod primes ≤ y₀) to remove the Murphy-α part. Both
reviewers showed this does not work as stated: (a) the residual D^{−2α}/y₀ ≤ D^{−(j−1)α} needs y₀ ≳ D^{(j−3)α}, not
D^{jα}; (b) the number of residue classes fixed is e^{(d+1)θ(y₀)}, and for j ≥ 5, (d + 1)y₀ ≥ D^{δ_j + (j−3)α_j} ≥ D (since
(j² − j − 3)/(3j + 2) ≥ 1), comparable to the log-size of the whole field space, so the conditioning pins the field
(the pointwise OD/AD regime); (c) archimedean fluctuations are not removed by local conditioning. Conditioning on a
low-dimensional statistic (a truncated Murphy-type sum plus size data) is the only candidate and is unexamined.
**So depth ≥ 4, and α → 1/3, remain open for this route with unknown feasibility;** the draft's "sub-L-scale at 1/3"
remark is withdrawn (it also ignored that the bounded-depth search cost grows with depth).

## 4. Fresh randomness (Q3)

| construction | what it refreshes | why it does not remove the joint problem | status |
|---|---|---|---|
| fresh field per level | rational children (integers) | algebraic children cannot be exported (Lemma FB, Cor EX); without them R2/R3 cap at 2/5 | proved (FB, EX) / model (R2, R3) |
| fresh field per child | — | an ideal of K₁ is not an ideal of an unrelated K₂ (Cor EX: k₁ = Q) | proved |
| fiber field per parent (F + G·L₀·H) | the d − 1 coefficients not pinned by P₀ | fresh only relative to the ≤ d − 1 retained parents; all their children close in that one field; other rows in it are ordinary one-field rows; private per-parent fields meet Lemma CC | CC proved; consequence model |
| multiplier w per child | one coordinate T = Nw | w is a function of the field; varying it per node changes the field (Prop NB) | proved |
| constrained resampling per component | d + 1 − (retained points) | same budget per algebraic component (nested §4) | proved algebra |
| Prop MF fields sharing g | algebraic side at one rational candidate | replaces the parent's algebraic side; does not export a child; rational values still share g | proved (MF), structural |
| union bound over children (§3) | nothing | needs per-child precision η ≍ 1/m instead of joint control | proved reduction (UB); PC₃ open |

**Lemma CC (singleton columns) [proved, trivial].** In the per-field square certificate class (escapeR2 §1), if a
factor-base prime ideal 𝔭 of a field K occurs with odd exponent in exactly one row of the collected matrix, that row
lies in no dependency. *Consequence (model, weakened after review).* A field used by only a few rows needs every
odd-exponent factor-base ideal of those rows to recur among them. Restricting to ideals of norm in (B^{1/2}, B] (of
which a typical smooth value has ≍ u_B, and which are not forced to coincide by small-prime bias), a row of a field
used by r rows has a private such ideal unless r is comparable to the number of active columns. So per-parent
private fields cannot be used as such; but a fiber field can of course serve further rows that are not retained,
which then get no fresh randomness. Neither this nor the row target M_rows ≥ C_col + J + 12 is a universal obstruction
to directed recombination (codex #10; proof critic minor 3, which also corrected the label count to ≈ π(B)).

**The rational side is the tightest budget** [proved; codex blind]: rank(G(P_i))_i ≤ 2, and 1 after conditioning on
G(P₀). So even three rational evaluations need a correlation estimate, whatever d is. In the 3/8 chain the rational
sides of all m level-2 points are B-smooth, so a joint estimate for complete certificates involves m + 1 ≈ d^{3/2}
rational forms in two variables, far more overdetermined than the algebraic block. Prop UB's per-child statements
involve three.

## 5. The barrier, sharpened (Q4)

**N&S for the pooled design class** [escapeR2 Theorem NR-q; codex blind (9)]: the class runs in L[α, O(1)] iff
#{q ∈ J : p_q ≥ e^{−cΛ_α}} ≥ (log N)^{−C₀}#J for some constants (T_q), equivalently an independent pool of e^{O(Λ_α)}
fields covers a polylog-dense root set with constant probability (sufficiency by amplification; necessity by
conditioning on a good pool and a union bound over its fields, codex). For the one-field design 2L the analogue is
γ₂ = Pr_{ℓ,A}[N₂ ≥ 1] ≥ (log N)^{−O(1)}, which is sufficient (Theorem 2L's proof) but **not necessary**: if a fraction
e^{−cΛ} of fields had δ_A = 1 and the rest δ_A = 0, finding a productive field costs L[α, c] while γ₂ = e^{−cΛ} (codex
#9). Necessity must be phrased through field-selection and per-field success probabilities. (The necessity direction
of T_q also rests on lv5unif's (R1), which is labelled heuristic.)

**Sufficient analytic inputs now known for 3/8.**
- LS_{3/2} (escapeR2): targetwise first and second moments of complete two-level certificates, ratio ≤ L^{O(1)}; ≈ d^{3/2}
  dependent algebraic forms in d variables and ≈ d^{3/2} rational forms in 2 variables; main term L[3/4]^{−Θ(1)}.
- PC₃ (this note): 3-point conditional moments, relative precision ≍ D^{−3/8}. Rational block = Hyp RB (ternary
  smooth equation with x^{o(1)} coefficients, family-averaged); algebraic block = BV-type 3-point lattice-coset
  asymptotics.
Neither implies the other as far as we can see. PC₃ has fixed geometric complexity (3 points, a complexity-1 rational
pattern) and polylog precision; LS_C has a growing number of dependent forms but only L^{O(1)} precision and needs no
conditional normalisation per incidence. Which is easier is a judgment, not a theorem.

**Necessary for the union-bound route** [proved: Lemma MP]: relative precision ≤ η for the child-box moments.
**Not known to follow from**: level-of-distribution hypotheses (EH for smooth numbers included). We know no derivation;
black-box AP statistics do not determine ternary additive counts for general sets, though the multiplicative
structure of smooth numbers (bilinear/dispersion methods) might in principle bridge the gap (proof critic minor 7). Codex blind gives a precise AP-statistics countermodel
for *fixed* searches (delete the ≤ W queried values: AP counts change by ≤ W, continuations vanish); for the averaged
PC₃ the obstruction is the standard one that AP statistics do not determine additive-equation counts. We know no
implication from GRH to RB: Lagarias–Soundararajan and Drappeau 2013 (GRH, coefficient 1, the latter covering
y = L[3/8] with relative error O(1/u)) do not address coefficients up to x^{o(1)}, conditional weights or covariance.
**Implied by**: the random model with local factors (heuristic, margin D^{3/8}, assuming uniformly bounded
perturbation variance); nested's AD (pointwise, much stronger); a covariance extension of Harper's theorem to
(i) coefficients up to x^{o(1)}, (ii) conditional weights from the algebraic block, (iii) relative covariance error
D^{−3/8−ε} ≍ ((log log x)/log x)^{1/2+ε} for the rational value scale x (the draft's "(log x)^{−3/8}" was too weak;
codex #5), together with a density-relative BV for the algebraic block.
**Implies**: nothing named that we could find (it is field-averaged; it says nothing about a fixed polynomial).

## 6. Experiments (pre-registered, `notes/multidescent_prereg.md`, commit `80d25ee`)

All toy: 48-bit balanced semiprimes (seeds 1–6), two main fields each (12 fields), d = 3 or 5. They test the shape of
the children's law and of conditional goodness; they support nothing rigorous and nothing asymptotic.

### 6.1 E1: the children's law (12 507 accepted level-1 points; 29 058 algebraic, 3 724 rational children)

Ratios (empirical share)/(reference share) in 4 equal-log-width sub-buckets of each window (2^j, 2^{j+1}]:

| window | alg. children | vs uniform labels | vs 1/q labels | rat. children | rat. vs 1/q |
|---|---|---|---|---|---|
| 2^11 | 5263 | 1.23 1.07 0.92 0.86 | 0.94 0.97 0.99 1.10 | 1076 | 0.97 1.07 1.03 0.93 |
| 2^12 | 5568 | 1.23 1.11 0.95 0.81 | 0.94 1.01 1.02 1.04 | 909 | 0.96 1.01 1.02 1.02 |
| 2^13 | 5815 | 1.18 1.07 0.97 0.86 | 0.90 0.97 1.05 1.10 | 849 | 1.11 1.04 0.98 0.86 |
| 2^14 | 5997 | 1.27 1.10 0.96 0.80 | 0.96 0.99 1.03 1.02 | 638 | 1.07 1.00 1.09 0.84 |
| 2^15 | 6415 | 1.28 1.07 0.92 0.84 | 0.97 0.97 0.99 1.07 | 252 | 1.76 1.06 0.74 0.41 |

- **P1a held**: every algebraic sub-bucket is within [0.896, 1.103] of the 1/q label law (band [0.85, 1.15]).
- **P1b held**: domination against uniform-on-labels, max 1.28, min 0.80 (bands ≤ 1.6, ≥ 0.6).
- **P1c failed** in the top window: rational children in (2^15, 2^16] are depleted towards the top (bucket counts
  114, 67, 46, 25; only the first has ≥ 100 children, and its ratio 1.76 exceeds 1.3). Explanation (post hoc, the
  §2.1 caveat): the rational cofactors G/ℓ are only ≈ 2^{17}–2^{20}, so a child q ≈ 2^{16} leaves a cofactor X/q of a
  few units and Ψ(X/q, y) collapses; the 1/q law omits the q^{1−β}Ψ(X/q)/Ψ(X) factor. The other windows hold
  (0.84–1.11). The algebraic values (≈ 2^{60}) are far from this edge.

### 6.2 E2: rational children in a pool of 24 fresh fields vs log-matched uniform primes (1 679 children)

| statistic | children | uniform (q′ ∈ [q/1.05, 1.05q]) | difference, 95% CI (parent clusters) | field-cluster CI (post hoc) |
|---|---|---|---|---|
| mean pool rate p_q | 0.1897 | 0.1868 | +0.0029 [−0.0023, +0.0082] | [−0.0019, +0.0076] |
| fraction p_q = 0 | 0.0208 | 0.0280 | −0.0072 [−0.0173, +0.0024] | [−0.0174, +0.0029] |

- **P2a held** (CI ⊂ [−0.03, 0.03]); **P2b held** (CI contains 0, ⊂ [−0.06, 0.06]). Matching on log q to ±5% removes
  the +0.04 offset nested E1 saw with ±25% windows (mean log₂q: 13.265 vs 13.267). Every parent had exactly one
  rational child at these parameters, so E2 does not test joint behaviour (none is needed by Prop RC).

### 6.3 E3: algebraic children under the conditional (fiber) law, and the rational block on lines

720 parents (d = 3, fiber rank 2; 81 fiber fields each) and 480 parents (d = 5, fiber rank 4; 81 of 625).

| statistic | d = 3 | d = 5 | band | verdict |
|---|---|---|---|---|
| P3a fiber-algebra violations | 0 | 0 | 0 | held |
| P3b within budget (m ≤ d − 1): R_joint (children) | 1.039 [1.03, 1.05] | 1.032 [1.021, 1.044] | — | — |
| P3b within budget: null R_joint | 1.036 | 1.012 | — | — |
| P3b alg/null, 95% CI | 1.003 [0.988, 1.017] | 1.020 [0.997, 1.041] | ⊂ [0.85, 1.2] | held |
| P3c overdetermined (m ≥ 3, d = 3; 144 parents): R_joint | 1.65 [1.45, 1.88] | (no parents) | — | — |
| P3c alg/null | **0.418** [0.283, 0.692] (null R_joint 3.94) | — | ∈ [0.67, 1.5] | **failed** |
| P3d rational block on lines: median excess rel. variance | −0.007 (1 334 boxes) | 0.052 (943) | [−0.1, 0.15] | held |
| P3e fiber block: median excess rel. variance | 0.123 (89) | 0.035 (250) | [−0.1, 0.4] / [−0.1, 0.2] | held |
| P3e Pr_H[N = 0]/e^{−μ}, median | 1.016 | 1.008 | [0.8, 1.25] | held |

- **P3c failed.** In the overdetermined stratum the children are *less* positively associated than the null, not
  more: R_joint 1.65 against the null's 3.94. Post hoc (`src/multidescent/posthoc.py`, not pre-registered): the toy
  fiber fields' heights grow with max|h_i| (the added term G·L₀·H has coefficients ≈ M·|P₀|·|h|), and height
  heterogeneity makes all targets succeed together (Jensen), the more so the more targets — more for rational
  targets, whose values F_H(P′) at P′ ∈ Λ_q(g) are larger. Computing J and Π *within* classes of equal max|h_i|:
  m ≤ 2 (564 parents): children 1.020, null 1.016, ratio 1.004; m ≥ 3 (138 parents): children 1.154, null 1.163,
  ratio **0.992**. So, once field height is held fixed, overdetermined children (m > fiber rank) are as close to
  conditionally independent as matched non-overdetermined targets. The pre-registered null was a poor control.
- **P3d/P3e.** Child-box counts are close to Poisson under both conditional laws: on the rational line g + jP₀ the
  excess relative variance is ≈ 0 (−0.007, 0.05), in the algebraic fiber 0.12 (d = 3) and 0.035 (d = 5), and the
  zero-probability matches e^{−μ} to 1–2%. At toy size the local factors are large (log B ≈ 8–9), so the excess is
  expected to be larger here than the asymptotic O(D^{−2α}); its decrease from d = 3 (B = 2^12) to d = 5 (B = 2^13)
  is consistent with that but two points are not a trend.
- **Reading (fixed in advance).** The toy shows nothing like a joint-control obstruction: children's goodness
  factorises up to field quality, also beyond the fiber rank, and the conditional counts that Prop UB needs to
  concentrate are nearly Poisson. This is a toy check of the premises of both Hyp LS_C and Hyp PC₃, not evidence for
  either asymptotically.

## 7. Literature (checked by the Opus literature critic against primary sources; all attributions confirmed or corrected)

- **Smooth solutions of a + b = c.** Harper, Compositio 152 (2016) 1121–1158, Cor. 1 (unconditional, y ≥ log^K x;
  minor arcs (log x)^{−4}, error O(log(u+1)/log y) from the simplified main term). Lagarias–Soundararajan, Proc. LMS
  104 (2012) 770–798 (GRH, (log x)^{8+δ} ≤ y ≤ e^{(log x)^{1/2−δ}}); also JTNB 2011 ("xyz conjecture"). Drappeau, MPCPS
  154 (2013) 439–463, Thm 1.3 (GRH, all y ≥ (log x)^{8+ε} with a mild condition, error O(1/u)); Drappeau, Canad. J. Math.
  67 (2015) (unconditional, y ≥ exp(c√(log x) log log x)); de la Bretèche–Granville, "Densité des friables" (Bull. SMF;
  y ≥ exp((log x)^{2/3+ε})). All have fixed or bounded coefficients; no covariance or growing-coefficient version found.
- **Linear patterns / binary forms.** Lachand, arXiv 1609.08872 (bounded u, bounded coefficients); Matthiesen–Wang,
  arXiv 2211.16892 (Ψ(N, y) > N/log^8 N, bounded coefficients; inapplicable at u ≍ D^{3/8}); Balog–Blomer–Dartyge–
  Tenenbaum, Comment. Math. Helv. 87 (2012) 639–667 (y = x^{g−2+ε}; x^{1/√e+ε} for cubics; order of magnitude).
- **BV for smooth numbers.** Granville 1993; Fouvry–Tenenbaum 1996; Harper, arXiv 1208.5992, Thm 1 (preprint);
  Drappeau, Compositio 151 (2015) 828–862 (Lemme 5; Théorème 1, level 3/5 − ε for residues a₁ā₂); Pascadi, arXiv
  2505.00653, Thm 1.5 (level 5/8 − ε, fixed residue). Hildebrand–Tenenbaum, Trans. AMS 296 (1986), Thm 3 (Ψ(cx, y) =
  Ψ(x, y)c^{α}(1 + O(1/u + log y/y)), 1 ≤ c ≤ y; used in Lemma CL, which needs z₂/q ≥ B₂).
- **Descent.** Granger–Kleinjung–Zumbrägel, Trans. AMS 370 (2018) 3129–3145 (Thm 3, trap removal); Kleinjung–Wesolowski,
  J. AMS 35 (2022) 581–624 ((pn)^{2 log₂ n + O(1)} for all p, n; hence L(α) for p = L(α), α < 1/2: rigorous multi-level
  descent exists in small and medium characteristic, not only fixed characteristic); Enge–Gaudry–Thomé (J. Math.
  Cryptol. 2011, heuristic), Commeine–Semaev 2006, Joux–Lercier 2003, Guillevic arXiv 1505.07553 — all heuristic, with
  implicitly independent children; none uses a per-child probabilistic union bound with a precision requirement.
  No rigorous multi-level descent for prime fields or Z/N found.
- **Factoring.** No rigorous or GRH-conditional factoring of arbitrary N below L[1/2] found; Lenstra–Pomerance, JAMS 5
  (1992), L[1/2, 1 + o(1)]. Lee–Venkatesan, JNT 187 (2018): L[1/3] congruences of squares, possibly trivial.
- This is a search result, not a novelty claim; "no prior formulation found by us" of the conditional 3-point
  reduction (Prop UB).

## 8. Reviews

All valid findings are applied in the text above.

**GPT-6.1-sol blind round `codex/multidescent_1`** (xhigh, launched before any analysis here). Verdict OPEN, no
α < 2/5. Independently: α_j = (j+1)/(3j+2); the exact children's law (1); the fixed-field domination obstruction (4);
the removal-sieve criterion (5) (it judged the needed small bad fraction unavailable — answered in §2.3); rank
(F, G) ≤ min(s, d+1) + min(s, 2), "even three rational evaluations require a correlation estimate"; T_q ⇔ pooled
coverage; an AP-statistics countermodel for fixed searches; Pascadi, Martin, BBDT, Harper, LS, Lachand, Matthiesen–Wang.
It did not find the union-bound reduction.

**GPT-6.1-sol adversarial round `codex/multidescent_2`** (xhigh, on v1): **MAJOR REVISION; no fatal flaw; Prop UB
sound.** Issues and adjudication:

| # | severity | finding | action |
|---|---|---|---|
| 1 | MAJOR | U′ needs a hybrid (B-smooth moduli, B₂-smooth values) restricted BV; the cited substitution fails both ways | accepted: U′ relabelled OPEN, requirements listed (§3.3) |
| 2 | MAJOR | box sizes: CS at B gives an L[1/2] box; S₂ unspecified | accepted: added to U′ |
| 3 | MINOR | "3 points" = 3 geometric points, ≤ 6 evaluations; §2.4 is a fibration with weights | accepted |
| 4 | MAJOR | algebraic BV reduction not established (density-relative, Harper's large-sieve term, Pascadi fixed residue); quantifier mismatch of pointwise PC₃ | accepted: reduction marked open; PC₃ restated in the weighted form codex proposed |
| 5 | MAJOR | precision D^{−3/8} is a typical-count prediction (cap gives D^{−1/2}); covariance O(D^{−2α}) needs uniform bounded variance; "(log x)^{−3/8}" wrong | accepted |
| 6 | MAJOR | Lemma CL's bad-moduli step invalid; RC sound conditionally and answers the blind objection | accepted: CL re-outlined via W-U (§2.2), marked UNVERIFIED |
| 7 | MAJOR | polylog per-child loss is affordable as work, not as root density (km loss) | accepted (§3.1) |
| 8 | MAJOR | local-data conditioning for depth ≥ 4 pins the field; 1/3 extrapolation unsupported | accepted: §3.5 rewritten, claim withdrawn |
| 9 | MAJOR | R3 not invalidated (certificates still 1 + m points); γ₂ not necessary (counterexample) | accepted (§3.4, §5) |
| 10 | MAJOR | Lemma CC's consequence unsupported (diffuse-label assumption; fiber fields can serve other rows) | accepted: consequence weakened (§4) |
| 11 | MINOR | Plücker coefficient ℓ·det(P′, P″); determinant sizes are upper bounds; Harper Cor. 1 error O(log(u+1)/log y) at exactly our scale; LS range excludes L[3/8]; GKZ trap removal is a union-type precedent; KW covers growing characteristic | accepted (§3.2, §3.4, §7) |

Codex's subjective probability that PC₃ for the actual design is provable by known methods in a focused effort:
≈ 15%. Ours: 5–10% (U′ and the conditional normalisation are additional hurdles).

**Opus proof critic (read-only, on v1).** Prop UB ACCEPT; Theorem 2L ACCEPT WITH FIXES (U′ not routine); Prop RC
ACCEPT WITH FIXES (constants checked); Lemma CL MAJOR REVISION (bad moduli in an upper bound; fix via W-U); Lemma CC
ACCEPT WITH FIXES (lemma trivial and correct, consequence's count wrong); model formulas ACCEPT; §3.5 MAJOR REVISION
(wrong cutoff y₀ ≳ D^{(j−3)α}; conditioning pins the field for j ≥ 5; archimedean fluctuations). Minor: primitive boxes
mod ±1 (perfectly correlated multiples), q ∥ F(P′) in Good, Chebyshev constant 3η/4, K ≥ 1 in RC, density 2ε′,
independence assumption in the per-child loss bookkeeping, necessity claims rest on (R1), EH wording, model.py
docstrings. All applied.

**Opus literature critic (read-only): OK WITH FIXES.** All attributions substantially correct. Applied: Harper's
O(log(u+1)/log y) is main-term bookkeeping, not a minor-arc limit (so the precision caveat is about bookkeeping, and
the real difficulties are coefficients, weights and averaging); Drappeau 2013 (GRH) covers y = L[3/8], correcting the
draft's "outside the GRH range"; LS's coefficient remark is one sentence; Drappeau 2015 Lemme 5 and Théorème 1 added;
GKZ precedent reworded as formal; §7 filled with the missing references.

## 9. Reproduction

```
.venv/bin/python -m pytest -q tests/test_multidescent.py
.venv/bin/python -m src.multidescent.experiments E1 E1   experiments/raw/multidescent_e1.jsonl 6
.venv/bin/python -m src.multidescent.experiments E2 E2   experiments/raw/multidescent_e2.jsonl 6
.venv/bin/python -m src.multidescent.experiments E3 E3d3 experiments/raw/multidescent_e3d3.jsonl 6
.venv/bin/python -m src.multidescent.experiments E3 E3d5 experiments/raw/multidescent_e3d5.jsonl 6
.venv/bin/python -m src.multidescent.analyze E1 experiments/raw/multidescent_e1.jsonl   (E2, E3 likewise)
```

## 10. Candidate card

- **Mechanism.** Two-level Dixon/NFS descent in one field (α = 3/8), children closed by a union bound with per-child
  concentration, instead of counting complete certificates.
- **Model.** Theorem U's design with d ≍ D^{1/4}, k ≍ D^{3/8}, child bound L[1/2], B = L[3/8].
- **Result.** Prop UB (proved); Theorem 2L: L[3/8] conditional on PC₃ and U′ (both open); Prop RC (rational children
  transfer, conditional on U_Q and Lemma CL); Lemma CC (trivial); α_j = (j+1)/(3j+2) (model).
- **Failure scope.** No unconditional α < 2/5. Fresh-randomness constructions fail for algebraic children (FB, EX,
  CC, NB). The high-point joint estimate (LS_C) is not needed by union-bound proofs, which need PC₃-level precision
  instead (Lemma MP); certificates remain large (R3, OD unaffected). Depth ≥ 4 is open for this route.
- **Next step.** Hyp U′ (hybrid restricted BV with B-smooth moduli and B₂-smooth values; mixed box scan), then Hyp RB:
  a covariance form of Harper's a + b = c with coefficients up to x^{o(1)} and conditional weights, beating the
  secondary-term scale D^{−3/8}; and the algebraic 3-point block by density-relative BV.
