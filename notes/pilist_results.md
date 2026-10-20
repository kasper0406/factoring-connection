# pilist: the product-listing primitive Π — where the L[1/3] limit bites

Agent `pilist`, phase 18, 2026-10-20/21. Revised after codex pilist_2 (MAJOR REVISION), the Opus proof critic
(MAJOR REVISION) and the Opus literature critic (§8).

**Files.**
- Code: `src/pilist/`.
  - `targets.py`: cost targets for Π (§1.3).
  - `decoder.py`: Prop JR bounds (unweighted and weighted), the bracketed count model, the HG/Boneh row
    decoder, and the row sieve harness.
  - `pencil.py`: split-polynomial enumeration, the line/pencil finder, and Lemma TL classes.
  - `tails.py`: the Prop TAIL first-moment bound and the F₂ tail simulator.
  - `run_experiment.py`: E1–E3.
- Tests: `tests/test_pilist.py`.
- Pre-registration: `notes/pilist_prereg.md`, commit c78c6bc, made before any registered run.
- Config: `experiments/configs/pilist.json`.
- Raw results: `experiments/raw/pilist_E{1,2,3}.jsonl`.
- Codex (GPT-6.1-sol): blind round `codex/pilist_1.{md,out}`; adversarial round `codex/pilist_2.{md,out}`.

**Labels.**
- PROVED.
- PROVED from KNOWN (a cited theorem is used as a black box).
- model-PROVED: an exact consequence of a named model.
- KNOWN: cited.
- HEURISTIC: the assumption is named.
- measured · open · assessment.

Novelty is never claimed; "no prior statement found by us" is the most we say.

**Notation.**
- x = log N, and F_α = x^α (log x)^{1−α}, so L[α, c] = exp((c + o(1))F_α).
- Factor base B = L[β]; box area U = E²; λ := log(C u log u)/log B, the per-unit-log cofactor gain of smoothenum
  Prop LQ step 1.
- For a decoder call: V is the size of the region it covers, P = ∏ q^{e_q} its modulus (factor-base primes),
  Λ = log P, S its threshold, and d(z) = gcd(f(z), P).

---

## 0. Verdict

**Plain language.**
- We found **no construction** and **no unconditional barrier**. An unconditional barrier would be a
  superpolynomial lower bound for an explicit arithmetic search problem, which current techniques cannot give.
- We did locate where the limit bites, in three conditional, method-level results.
  1. **Certificate decoders stop at the Johnson radius.** Take any Coppersmith / Howgrave-Graham / Boneh /
     Guruswami–Sahai–Sudan CRT list decoder: any dimension, any modulus built from factor-base prime powers, any
     weights. Within its guaranteed (Johnson-type) range, every certified output has a known divisor at least as
     large as the region searched.
     - Batching all those divisor lattices into one call is then worth only L^{o(1)}. This holds in the
       random-divisibility model on oblivious regions, for primitive points, in NFS scaling (the congruence
       instances are excluded) (Prop JR, model-PROVED).
     - So Boneh's L^{o(1)} is a barrier for that whole class, not an artifact of how Boneh used it.
     - Doing better needs output-sensitive listing beyond the Johnson radius. That would *be* a partial lister,
       i.e. the open problem itself.
  2. **Adaptive Θ(u)-tail elimination saves at most an L^{o(1)} factor in relations.** With boundedly many forced
     tail primes per relation, a single dependency still needs ≥ 2^{−8}/(kp) = L[β]^{1−o(1)} relations (Prop TAIL,
     under the strengthened randomness hypothesis H_tail⁺). What could still help:
     - sub-quadratic sparse linear algebra;
     - H_sel-type selection;
     - forcing by other means;
     - failure of H_tail⁺ itself.
  3. **For Π_β itself the target is extreme.** The congruence instances must cut the per-output cost from L[1 − β]
     (box sieve) to L[1/3, O(1)] even just to beat MNFS. Only the NFS-shaped listers have fixed-fraction
     thresholds.
- **smoothenum's crux H_sel is not equivalent to the primitive** (a draft claim, withdrawn after review). Its gain
  form H_sel^gain is exactly the hypothesis LQ uses. A violation of H_sel^gain with fresh primitive outputs, i.e. of
  H_sel by ≥ (M/U)·L[θ′r]^{−1+ε}, is a sufficient criterion for a speedup.
- **PTE O2.**
  - The naive three-term rigidity is false: explicit identities at M = 3 (exact) and M = 4 (constant remainder).
  - Even a true Lemma R-shaped version could not bound O2, because of Behrend-type parameter sets.
  - The right statement is a count. A nondegenerate pencil has ≤ (16/3)(cB^{1/2}(log B)^κ + 1) + min(2e − 2, 2B + 1)
    split members, height-free (Lemma PEN; codex proved it independently).
  - A pair of presentations with remainders of degree < k has ≤ C(2B+k, k)·(that) nondegenerate members (Prop
    O2-count). Degenerate W-pencils can add O(B²).
  - That is bounded amplification, not exclusion.

**Main results and labels.**

| Result | Status |
|---|---|
| Π canonical: split values of binary forms on a box; Table 1.1 targets | PROVED containments; model-PROVED / numerical costs |
| Remark SEL (H_sel vs partial listing; one-primitive-point lemma) | PROVED, elementary; the equivalence is withdrawn |
| Lemma CS (certified ⇒ d ≥ V) | PROVED |
| Prop JR (weighted and unweighted Johnson-radius certificates, Γ_cert = L^{o(1)}) | model-PROVED: H_div on oblivious regions, primitive points, NFS scaling. Proof: codex, checked by us. Not a proof of H_sel |
| g_n(μ) ≤ μ² (Herrmann–May ranges lie inside Johnson) | PROVED (AM–GM) |
| Prop TAIL (R ≥ 2^{−8}/(kp) for one dependency) | model-PROVED under H_tail⁺ and LQ |
| Lemma TL, Lemma PEN, Prop O2-count, the TR counterexamples, the Behrend limitation | PROVED / PROVED from KNOWN (BCN; Behrend) |

**Experiments (pre-registered, commit c78c6bc; every failure listed).**
- **E1.** P1, P2 and P5 held.
  - **P3 failed** at M = 4: no exact coprime small-c pencil with e = M.
  - **P4 failed as scored**: 2 flags, both G = 0 triples with moving degree 2. A disclosed post-hoc rescoring
    removes them.
  - Exact nondegenerate pencils carry ≤ 3 members up to M = 10.
- **E2** (c60 rational side). P6 (0 misses in 272 HG calls), P8 (≤ 0.77 certified outputs per call) and P9 held.
  - **P7 failed** in 1 of 8 cells (×2.37). Post hoc: non-primitive multiples of (30, 1). The primitive density
    1/(q + 1) fits.
  - Pooled Γ̂ ≈ 0.30: a call is worth less than one box point.
- **E3.** P10 and P11 held: legal adaptive strategies have median R*/n ≥ 0.917. The illegal oracle reaches
  0.79–0.95.
- One crash-and-fix deviation in E1 is disclosed.

**Assessment (subjective).** Probability that Π_β is solvable for some β < 1/3: below 1% (ours); about 5% (codex,
blind). The results are conditional and method-level; they do not calibrate this number.

**Where to push next.** An output-sensitive lister beyond the Johnson radius for the factor-base CRT code, i.e. for
the zero word f(t) mod q. Plain-RS capacity decoding (Brakensiek et al. 2026, as reported) shows that
beyond-Johnson decoding without auxiliary data is possible for RS codes. Whether anything like it exists for CRT
words, at vanishing rate and with output-sensitive huge lists, is the sharpest open form of the primitive.

---

## 1. The canonical problem Π and its cost targets (Task 1)

### 1.1 Definition

**Π (split values of binary forms on a box).**
- *Input.* N, and data computable from N:
  - binary forms f₁, …, f_r ∈ Z[a, b] of degrees δ_i;
  - multiplicative bases 𝔛_i ⊂ Z, which are sets of integers: primes ≤ L[β], or shifted pieces such as m − a;
  - a box Ω ⊂ Z².
- *A witness* is a point z ∈ Ω such that every f_i(z) is 𝔛_i-split (a product of elements of 𝔛_i), together with
  the factorisations.
- *Task.* List R ≥ K^{1+o(1)} witnesses (K = Σ|𝔛_i|) at amortised cost τ per witness. Their parity vectors must
  have F₂-rank ≥ (1 − o(1))·#columns. Charge everything: building the forms and bases, setup, extraction, duplicate
  removal.
- An optional *window side* j allows |f_j(z)| ≤ V (bounded by construction) in place of f_j being split. The
  window side must then be smooth for the relation to be usable; this is the one-sided case.

**The three primitives are instances [PROVED, elementary].**

| Primitive | Forms and box | Bases |
|---|---|---|
| Π-nonexact (one-sided, View A at m = 0) | f₁(b, v) = bN + v on [1, U/N] × [−V, V]; side 2 is the window v | primes ≤ L[β] |
| Π-nonexact (View A, pieces m − a) / Π-split (ptecons item 3) | the same form bN + v | the shifted base {m − a : \|a\| ≤ B}. A split product ∏(m − a_l) is the lift P_A(m) of a split polynomial |
| Π-nonexact (doubly split) | the same form, with v also split | primes, or m − a |
| Π-NFS rational / algebraic / doubly | f₁ = a − bm (rational), f₂ = F(a, b) (homogeneous, degree δ) on [−E, E] × [1, E]; one side split, the other tested; or both split | primes / prime ideals |

- *Why the congruence form is a binary linear form.* A product u = ∏x_a ≡ v (mod N), with |v| ≤ V and |u| ≤ U,
  is the same thing as a point (b, v) = ((u − v)/N, v) of the box at which bN + v is split. Conversely, every such
  point gives a product.
  - This is the bijection of smoothenum Prop GEN, here at the modulus N.
  - The rational NFS side is the same statement at the modulus m ≈ N^{1/(d+1)}.
- *View B.* Residues μ − a ∈ Z/N have no integer lift of size < N, so View B is the form bN + v with the base
  {μ̃ − a} of lifts. Its pieces are of size N; only the doubly-split version is meaningful.
- *Direction of the containments.*
  - Each primitive is an instance of Π: the forms, bases and box above define exactly its witness set.
  - The algebraic NFS side needs prime-ideal bookkeeping (norms factor; ideals are the columns). That is extra
    data in the factorisation, not a different listing problem.
  - No reduction in the other direction is claimed.

**Scope.**
- Π is a listing problem only. Factoring also needs fruitful square dependencies, linear algebra and (one-sided)
  the smoothness of the window side (nonexact §6; heurnfs H; Croot–Granville–Pemantle–Tetali 2012).
- Literal listing costs ≥ 1 per output (codex). "Below one unit per output" means, here and in the task, below the
  box sieve's cost per output.

### 1.2 Two regimes, one structural fact

For the congruence instances the box is huge.
- To hold L[β] witnesses, bN + v must be split; it has size ≥ N, so its smoothness probability is
  L[1 − β]^{−Θ(1)}. Hence the box has area ≥ L[1 − β].
- So the box sieve costs **L[1 − β] per output**, while Π_β asks for **L[β, O(1)] per output**.

For the NFS instances the box has area L[1/3, 2e] and the sieve costs L[1/3, r] per rational-smooth output.

**Consequence [model-PROVED, heuristic H].**
- On the NFS instances, a fixed fractional saving of the box sieve's per-output cost beats MNFS: any θ > 0 on the
  rational side, and θ > θ* ≈ 0.03 / 0.022 on the algebraic / doubly side (numerical).
- On the congruence instances (Π_β proper), every target, *even only beating MNFS*, needs the per-output cost to
  drop from L[1 − β] to L[1/3, O(1)]. In the partial-listing language that is θ = 1 − O(F_{1/3}/F_{1−β}) → 1, i.e.
  essentially full output sensitivity.
- So, as a target and against the box-sieve baseline, Π_β demands far more than the NFS-shaped partial listers.
  This is a comparison of targets, not a proved complexity separation. Its natural home is the doubly-smooth lister
  of smoothenum T7/T8.
- The factoring consequence of either target also assumes random window smoothness and useful rank (codex #12).

### 1.3 Cost targets

**(A) Congruence instances, level α** (`targets.py`).
- Parameters: factor base L[α, κ]; window |v| ≤ L[2α, ν], i.e. log V = νF_{2α}; lister cost L[α, τ] per output.
- One-sided: v is smooth with probability L[α, −αν/κ]. The total exponent is
  - max(κ + τ + αν/κ, 2κ) one-sided;
  - max(κ + τ, 2κ) doubly split.
- Sanity check: α = 1/2, ν = 1, τ = 0 is Dixon's L[1/2, √2] (tested).
- The largest admissible τ for a target constant c [PROVED]:
  - τ_max = c − 2√(αν) if √(αν) ≤ c/2;
  - c/2 − 2αν/c otherwise;
  - doubly split: c (approached as κ → 0).

**(B) NFS instances (MNFS program of heurnfs/smoothenum).**
- A partial lister enumerates its side's smooth pairs at cost E²·V·L[−θ·penalty].
- The penalties: r (rational), a (algebraic), r + a (doubly).
- c_rat(θ) uses the closed form of smoothenum/codex. c_alg and c_both are 5-variable Nelder–Mead programs with e ≤ 3
  *enforced as a constraint*, with a running minimum against MNFS. They are numerical upper bounds, not certified
  optima. The draft's cap was not enforced (codex #11); the values below are the corrected ones. The algebraic and doubly listers enumerate per field (cost includes V);
  rational outputs are tested in all V fields.

**Table 1.1: what Π must cost, per output.**

| Target | Congruence Π (Π_β, Π-split, Π-nonexact) | NFS-shaped Π |
|---|---|---|
| beat MNFS 1.9019 at α = 1/3 | one-sided: τ < 1.9019 − 2√(ν/3). Window ν = 1 needs per-output cost < L[1/3, 0.747] (ν = 0.5: 1.085; ν = 2: 0.269). Doubly split: any τ < 1.9019. Box sieve: L[2/3, Θ(1)] per output | rational: any θ > 0 (c = 1.8998, 1.8912, 1.8804, 1.8467, 1.786, 1.6386 at θ = .01, .05, .1, .25, .5, 1; slope −0.21). Algebraic: θ > θ*_alg (c = 1.8867, 1.8497, 1.7347, 1.5247, 0.8246 at θ = .05, .1, .25, .5, 1; 1.9019 at θ = .01). Doubly: θ > θ*_both (1.8745, 1.8253, 1.6711, 1.3868, 0.5732) |
| any α < 1/3 | τ = O(1) at level α with ν = O(1): per-output L[α, O(1)] against the sieve's L[1 − α] | doubly listing with 1 − θ ≲ F_α/F_{1−2α} → 0 (codex eqs. 8–10), together with the doubly-smooth side |
| α = 1/4 | L[1/4, c] iff τ + √ν < c (when √(ν/4) ≤ c/2). E.g. ν = 1/4 gives τ < c − 1/2 | 1 − θ ≲ (1/8)x^{−1/4}(log x)^{1/4} (codex eq. 10; θ ≥ 0.9923 at x = 10⁶) |

- θ*_alg = 0.0303 and θ*_both = 0.0221 (bisection, 30 restarts, cap enforced). Codex's reduced program gives 0.0291
  and 0.0219. Below them, per-field enumeration costs more than MNFS's shared
  rational sieve.
- 0.8246 and 0.5732 are the e ≤ 3 values (codex: 0.8241 and 0.5716). The infima, attained only as e → ∞, are 0.7631 and 0 (smoothenum T6, T7).

---

## 2. Barrier attempts (Task 2)

### 2.1 What the crux is (Remark SEL; revised after both reviews)

**Lemma (at most one primitive point per large lattice) [PROVED, elementary].**
- For a rational-side lattice Λ_{M,R} with M > 2E², all box points are integer multiples of one vector: two
  independent box vectors would span a parallelogram of area ≤ 2E² < M = det.
- So Λ_{M,R} contains at most one *primitive* box pair up to sign.
- Caution (codex, critic): Λ_M can be nonempty with **no** primitive box point. Example: E = 2, M = 12, R = 5,
  where (−2, 2) ∈ Λ but (−1, 1) ∉ Λ. "Nonempty" must mean "contains a primitive point" throughout.

**Remark SEL [PROVED, elementary; the draft's "Prop EQ-SEL" equivalence is WITHDRAWN].**
- *Two hypotheses.* smoothenum's LQ theorem uses
  - H_sel: E[#new box points of Λ_M | 𝓕] ≤ C·U/M for M > U;
  - H_sel^gain: for every selector, (π(M)/T)·g(M) ≤ L^{o(1)}, where π(M) is the expected number of fresh primitive
    box points per selection and g(M) = ρ(u − log M/log B)/ρ(u).
  - H_sel implies H_sel^gain, but not conversely.
- (a) H_sel^gain is exactly the hypothesis LQ's proof needs. A partial lister with θ > 0 that acts through LQ
  operations must violate it. This is LQ's contrapositive, restated.
- (b) The draft's converse ("no partial lister ⇒ H_sel") is false. If revealed points may be re-selected, the
  ordinary θ = 0 box sieve already violates H_sel: it finds a smooth z at cost L[r] and then sets up Λ_{|G(z)|}.
  If they may not (smoothenum counts only new points), a listed point supplies nothing new. Either way the
  argument cannot separate partial listing from sieving (codex #6, critic M3).
- (c) Sufficient criterion. A selector that violates H_sel by φ(M) at moduli M = y^{θ'}, *with fresh, distinct
  primitive outputs*, gains φ(U/M)g(M) per unit cost. If φ ≥ (M/U)·L[θ'r]^{−1+ε}, that is an improvement with
  θ = εθ' (Dickman-gain assumption). It is a criterion, not an equivalence.

**What the natural models give** (codex blind round, checked):
- *Hidden fresh residues*: H_sel holds with C = 1 [model-PROVED]. But real roots are reductions of one public m,
  so the model hides the arithmetic the question is about.
- *Public unrestricted m*: a uniform H_sel is false (m = 1 + M with smooth M ≫ U puts (1, 1) in Λ_M). H_sel can
  only be an averaged statement.
- *Generic CRT/lattice manipulation with public residues*: a lower bound would be a superpolynomial lower bound
  for an explicit arithmetic search problem. We know no technique for that (assessment).

**Conclusion.** The crux can be attacked only in method classes, or as a conditional statement. §2.2 does both for
the largest method class we found.

### 2.2 Prop JR: Johnson-radius certificates batch lattices but gain only L^{o(1)} [model-PROVED; conditional]

**Scope (restricted after review).**
- NFS scaling: the call's region lies inside the NFS box, so V ≤ U with log U = O(F_{1/3}) = O(log B).
- The congruence instances are **not** covered. There log V/log B → ∞, and a single call with P ≈ d ≈ V
  certifies one point with gain g(d) = L[1/3, Θ(1)] (critic M2). That is harmless for Table 1.1, because the
  per-output cost stays L[1 − β].

**Certified outputs.** The class's output is *defined* as Cert = {z ∈ R : W(z) ≥ T}. Post-filtering is free;
uncertified returns are discarded, and the theorem says nothing about them. Uncertified returns are chosen
arithmetically (roots of h), so one could sit just below S with gain ≈ S^λ.
The draft let a call return an arbitrary superset; codex showed that permits Γ = V (codex #1, critic M5).

**Decoder class 𝒟_J (weighted; unweighted is w ≡ 1).** Guruswami–Sahai–Sudan's weighted agreement, with
nonnegative weights w ≢ 0 on the primes of P = ∏q^{e_q} (prime powers allowed):
- W(z) = Σ_q J_q(z) w_q log q, where J_q = min(v_q(f(z)), e_q);
- certified iff W ≥ T := √(log V · Q), with Q = Σ e_q w_q² log q.
- The unweighted case w ≡ 1 is Boneh/GRS/HG: log V ≤ (log S)²/log P.
- GSS's actual theorem (FOCS 2000, Thm 4) has an extra tolerance √(… + ε z_max²) and guarantees inclusion, not
  "iff". Its moduli are pairwise coprime, with Guruswami–Sudan multiplicities.
- Our prime-power score J_q = min(v_q, e_q) is therefore a *GSS-type extrapolation*, not literally a GSS
  coordinate. Unweighted, it is covered by Boneh/HG in gcd form.
- Lattice dimension, shifts and multiplicities are arbitrary; only the guaranteed range enters.

**Model H_div (stated after review).**
- *Regions.* Oblivious regions: rows or sub-boxes fixed before the call. Adaptive sequences of calls are allowed,
  with the hypotheses conditional on the transcript.
- *Points.* Primitive points only.
- *Divisibility.* For z uniform in R and q | P, the J_q are independent with Pr[J_q ≥ j] ≤ q^{−j}. The primitive
  density 1/(q^{j−1}(q + 1)) satisfies this (E2, §5.2).
- *Cofactor.* Given d(z) = d, the cofactor is smooth with probability ≤ P(y)·g(d), with g(d) ≤ d^λ (Dickman step).
- *What H_div excludes.* Arithmetic choice of the region. R = Λ_M ∩ box would violate it, and is H_sel itself
  (critic M1, codex #4). So **JR does not prove H_sel**. It proves that a Johnson-radius call, batching all
  divisor classes of P at once, cannot escape the H_sel-type count, *given* random divisibility on oblivious
  regions.

**Lemma S≥V and Lemma CS [PROVED].**
- Every certified output has d ≥ V. Proof: by Cauchy–Schwarz, W² ≤ (Σ_agree w² log q)·log d ≤ Q·log d, and
  W² ≥ T² = Q log V.
- In the unweighted case, S ≥ V^D follows directly from (log S)² ≥ D log V·log P ≥ D log V·log S.
- Up to a factor of 2 between V and Coppersmith's radius X = (V − 1)/2.

**Prop JR.** Assume NFS scaling, F := F_{1/3}, t = log V ≤ log U = O(F), λ ≤ 1/10 (true asymptotically), and
H_div. Then for every factor-base modulus P (any prime powers) and every nonzero weight vector,

  Γ_cert ≤ V·E[d^λ·1_{W≥T}] ≤ exp( λt + O(F/log F + F^{0.9}log F + u + log log B) ) = L^{o(1)}.

**Proof** (codex pilist_2's weighted proof, checked step by step by us; it replaces the draft's three-regime
proof, which the critic had found correct up to constants and which is kept in `decoder.py` as the numeric
Rankin bound).
1. Rescale w so that Q = t; the certified event becomes W ≥ t.
2. With s = 2(1 − λ), exponential Markov and independence give
   log Γ_cert ≤ t − st + Σ_q m_q(s w_q + λ), with m_q(a) = log E[q^{aJ_q}].
   Monotone domination by truncated geometric variables is legitimate because all exponents are ≥ 0.
3. *Heavy primes* H = {w_q > 1/4}. From Q = t, log P_H ≤ 16t. Also m_q(a) ≤ e_q(a − 1)₊ log q + log(e_q + 1),
   and (a − 1)₊ = (1 − λ)(2w − 1)₊ ≤ (1 − λ)w². So
   Σ_H m_q ≤ (1 − λ)Q_H + log τ(P_H) ≤ (1 − λ)t + O(F/log F).
4. *Light primes* (w_q ≤ 1/4): a ≤ 0.55.
   - Primes q ≤ Z = F^{1/4} contribute O(Z) uniformly in e_q.
   - For q > Z, Taylor-expand m_q around λ, using m_q′(λ) ≪ (log q)q^{λ−1} and m_q″(a) ≪ (log q)²q^{−0.4} for
     a ≤ 0.6:
     - base term Σ_{q≤B} m_q(λ) = O(u + log log B);
     - linear term O(√t), by Cauchy–Schwarz against Q = t;
     - quadratic term O(t(log Z)Z^{−0.4}) = O(F^{0.9} log F).
5. Sum: t − 2(1 − λ)t + (1 − λ)t = λt, plus O(F/log F) (divisors), O(F^{1/4}) (small light primes), O(√t)
   (linear), O(F^{0.9} log F) (quadratic) and O(u + log log B) (base). λt = O(λ log U) = O(log u). Every error term
   is O(u) or smaller (proof critic). ∎

**Reading.**
- Lemma CS says every certified output has a known divisor at least the size of the region. Each divisor class
  then holds ≤ V/d ≤ 1 expected points, the H_sel-type count.
- The decoder's only freedom is to batch the τ(P) classes into one call. The proof shows batching is worth
  L^{o(1)} for every weight profile: heavy primes are capped by Q = t, and light primes contribute nothing beyond
  their mean.
- The unweighted Johnson-radius class (Boneh/GRS/HG) is the special case w ≡ 1.

**Numbers (`src/pilist/sweeps.py`, the convex Rankin bound min_s V S^{−s}E[d^{s+λ}]).**
- Unweighted, z ≤ 10240 and powers ≤ 6: max log Γ = 3.13/4.87 (row/box) at 512 bits, 3.18/5.42 at 1024,
  3.98/4.07 at 4096, 4.53/−7.49 at 65536. Row maxima rise slowly (3.1 → 4.5). That is consistent with the proved
  exp(O(u + …)), which is L^{o(1)} but grows faster in u than smoothenum's (C u log u)^{v_box}.
- The maximiser is always log S ≈ log P ≈ log V: the single-lattice peak.
- Weighted, four families: max 5.39 at 1024 bits, < 0 at 65536.
- At E2 scale λ = 0.132 > 1/10, so E2's "JR bound" column is the numeric Rankin bound, not the theorem.

### 2.3 Is Boneh's L^{o(1)} a barrier or an artifact?

**A barrier, conditionally, for the Johnson-radius certificate class** (unweighted and weighted, any modulus,
prime powers, any dimension; model H_div; NFS scaling).
- It is not an artifact of Boneh's row-by-row use, of the choice of P, or of the weights.
- *Optimality of the method.* Chinburg–Hemenway–Heninger–Scherr (ASIACRYPT 2016) prove, by capacity theory, that
  Coppersmith's univariate full-modulus bound is optimal for auxiliary-polynomial certificates. This supports
  (J) as the method's real limit in the univariate case. It does **not** cover the divisor version or
  multivariate certificates (codex #15, critic M4). For those, (J) is the known guarantee, not a proved ceiling.

**An artifact relative to beyond-Johnson listing.**
- A useful decoder must list every z in a region of size V ≤ U with B-smooth part ≥ y^{θ'}, using the whole
  factor base (log P ≍ B). That is agreement fraction log S/log P ≍ x^{2/3}/L[1/3]: far below Johnson, above
  capacity. Doing it output-sensitively is a partial lister with θ = θ' (Remark SEL (c)).
- *Literature status (corrected).* The draft's "every capacity-approaching decoder needs folded or derivative
  structure" is withdrawn.
  - Brakensiek–Chen–Putterman–Zhang–Zheng, "Algorithmic list decoding of Reed–Solomon codes up to capacity"
    (arXiv 2609.08005, Sept. 2026; preprint, not yet refereed; verified by the literature critic). It list-decodes
    *plain* RS codes over prime fields up to radius 1 − ε, at rates up to (1 − θ)ε, extended to all constant rates.
    The decoder gets no derivative information; Hasse derivatives are hidden formal variables.
  - Companion: Jeronimo, arXiv 2609.05870.
  - It needs a prime field. Lists have size n^{O_ε(1)} and time is q^{ε^{−O(1/θ)}}. As the agreement fraction
    → 0, here ≍ x^{2/3}/L[1/3], that exponent blows up superpolynomially. It is not a CRT decoder, so it does not
    transfer to our vanishing-agreement, huge-list, output-sensitive regime.
- *The known integer beyond-Johnson mechanism needs extra structure.* Cohn–Heninger, "Approximate common divisors
  via lattices" (ANTS X, 2013), the integer analogue of Parvaresh–Vardy, goes beyond the Johnson-type bound only for
  several correlated inputs sharing one divisor. Our setting has a single word f(z) mod P.
- No CRT or number-field extension of the hidden-derivative method was found (literature critic).
- Cohn–Heninger AMC 2015 Thm 1.3 is the number-field β²/d Johnson-type bound, not capacity (draft misclassified
  it; per codex).
- *Precursors.* Boneh 2002 §4: the CRT list-decoding bounds "are insufficient for improving the basic quadratic
  sieve". Boneh's Thm 3.2 already bounds the number of N ∈ I with gcd(N, S) > T, the per-class count used in JR.
  Guruswami's thesis §12.3.2 makes the same remark. No published L-power impossibility for the Johnson-radius
  class was found by us or the literature critic.
- *Chinburg et al.* say explicitly that the unknown-divisor version "appears to pertain to joint capacities of
  many adelic sets, a topic which has not been developed". The divisor/CRT optimality question is open.
- Bundling primes into blocks yields (n/w)(u/n)^w → 0 agreeing blocks (codex blind, eq. 17).

**Comparisons.**
- *RS decoding hardness* [KNOWN]:
  - Guruswami–Vardy (IEEE TIT 51, 2005): maximum-likelihood decoding is NP-hard;
  - Cheng–Wan (FOCS 2004 / SICOMP 37 (2007); IEEE TIT 56 (2010), bounded-distance at (2/3 + ε)·distance): links
    to discrete logarithm;
  - Gandikota–Ghazi–Grigorescu (FOCS 2016): NP-hardness near n − k − O(log n/log log n).
  - Our direction is the same: beyond-Johnson output-sensitive listing for the factor-base CRT code ⇒ an L-power
    NFS gain. This is a consistency check. None of these results is average-case for the arithmetic words here.
- *Enumeration lower bounds* [KNOWN]:
  - Aono–Nguyen–Seito–Shikata (CRYPTO 2018): pruned enumeration;
  - Hanrot–Stehlé, "Worst-case Hermite–Korkine–Zolotarev reduced lattice bases" (arXiv 0801.3331): tightness of
    Kannan's n^{n/(2e)}. The draft wrongly cited CRYPTO 2007.
  - These are method lower bounds, like nonexact Prop LIN. They say nothing about selection among lattices, which
    is where JR operates.

---

## 3. Construction attempts (Task 3)

### 3.1 Whole-box multivariate decoders

**Linear forms (rational side).**
- Herrmann–May (ASIACRYPT 2008, Thm 4; under their Assumption 1, with ε slack): ∏X_i < P^{g_n(μ)}, where
  g_n(μ) = (n+1)μ − n + n(1 − μ)^{1+1/n}.
- **g_n(μ) ≤ μ² [PROVED].** With a = 1 − μ, μ² − g_n(μ) = a[(n − 1) + a − n·a^{1/n}] ≥ 0 by AM–GM (codex; checked
  to 10⁻¹⁴ on a grid).
- So the Herrmann–May decoder, and any 2-D decoder whose range is at most theirs, satisfies (J) with V = the box.
  Prop JR then applies: **Γ_cert = L^{o(1)}** [model-PROVED, conditional on HM's heuristic for the range].
- Caveats (critic M4):
  - HM is stated for inhomogeneous linear f. For the homogeneous G = Y₁a + Y₀b the algebraic-independence step
    degenerates (only a/b is determined) and needs a separate check.
  - Nothing here bounds *all* 2-D certificates.
- Not covered: whole-box decoders whose range exceeds HM's. None is known to us.

**Algebraic side (degree δ) and doubly-smooth decoders.**
- (J) has D = δ, a smaller range. Jochemsz–May (Assumption 1) treats known moduli, so its use for unknown divisors
  is an extrapolation of their strategy.
- The JR proof applies with ρ_q ≤ δ roots per prime, provided root counts behave on average as the prime ideal
  theorem predicts. The worst case ω(P) log δ is not o(F).
- HEURISTIC.

**Cost in smoothenum's model.** 1.9019 (gain L^{o(1)}). Codex's batch condition, ω − κ ≥ θr, fails for the class.

### 3.2 GS / PV / folded / multiplicity CRT decoding

- **Within Johnson** (GRS; Boneh; GSS with weights and multiplicities): Prop JR, L^{o(1)}.
- **Beyond Johnson.**
  - What would suffice is output-sensitive listing of {z : smooth part ≥ y^{θ'}} with the full factor base. By
    Remark SEL (c) that is a partial lister. So this open item is not a separate lever: within Johnson it is
    closed conditionally (JR); beyond it, it is the crux.
  - Known beyond-Johnson decoders (folded RS; multiplicity codes; plain RS, Brakensiek et al. 2026 preprint, §2.3) work at
    fixed rate, with n^{O_ε(1)} lists and time that blows up as agreement → 0. None is a CRT decoder for the vanishing-rate regime here.
- **Cost:** within Johnson 1.9019; beyond, HYPOTHETICAL c(θ') from Table 1.1.

### 3.3 Adaptive elimination with Θ(u) tail primes (Prop TAIL; restated after review)

**Model H_tail⁺ (strengthened after codex #8).**
- Band columns are the primes in (B^{1−ε}, B], with ε = C/log u.
- A B-smooth value of size B^u has (1 − e^{−C} + o(1))·u band factors, counted with multiplicity (Dickman saddle
  point; HEURISTIC). The model assumes **exactly k** random band entries per relation; the bound is for row
  length exactly k, and odd unmatched configurations contribute 0.
- Each relation has ≤ j forced band entries chosen from the history. Lattice index ≤ U gives
  j ≤ log U/((1 − ε) log B) = O(1).
- Entry by entry: given the entire exposure history, including earlier entries of the same relation, each random
  entry equals any fixed prime with probability ≤ p := C′u log u/B^{1−ε} = L[β]^{−1+o(1)}.
- Without the entry-wise condition the bound fails: rows (Q, Q, Q′, Q′) always cancel (codex).
- *Relations with fewer than k band primes* (proof critic) cost only an extra e^{−Θ(u)} = L^{o(1)}, so they are
  not excluded. Discarding them, or keeping only relations whose band count is ≥ k, amounts to a factor base
  B^{1−ε} = L[β]^{1−o(1)} on the remaining columns. The conclusion R ≥ L[β]^{1−o(1)} is unchanged.

**Prop TAIL [model-PROVED].** Under H_tail⁺, for any adaptive strategy (discards charged), and with
k ≥ 2j + 2:

  Pr[some nonempty subset of the first R relations cancels on the band columns] ≤ Σ_{s=1}^{R} C(R, s)·Σ_{w≤js} C(ks, w)·(ks − w − 1)!!·p^{(ks−w)/2}

and this is < 1/2 whenever **R ≤ 2^{−8}/(kp)**.

**Proof.**
- If S cancels, the XOR of its random parts equals the XOR of its forced parts, of weight ≤ js. So ≥ (k − j)s/2
  disjoint position pairs carry equal primes.
- For a fixed S and pairing, the probability is ≤ p^{#pairs}: condition at the later entry of each pair.
- Union bound, with C(R, s) ≤ (eR/s)^s, C(ks, w) ≤ 2^{ks}, (n − 1)!! ≤ n^{n/2} and ksp ≤ 1:
  term_s ≤ (js + 1)·[e(R/s)·2^k·(ksp)^{(k−j)/2}]^s.
- With R ≤ c/(kp), ksp ≤ cs/R, and since (k − j)/2 ≥ 1, the bracket is ≤ e·2^k·c^{(k−j)/2}. Using
  k − j ≥ k/2 + 1, this is ≤ e√c·(2c^{1/4})^k.
- At c = 2^{−8} it is ≤ e·2^{−4−k}, and Σ_s (js + 1)(e2^{−4−k})^s < 1/2 for j ≤ k/2. ∎ (critic m1)

**Constants.** c₀ = 2^{−8} is proved. The exact first-moment threshold gives R·k·p at the 1/2-level of:
- 2.48, 2.70, 2.72, 2.72 for j = 0, k = 4, 8, 16, 32;
- 0.64 → 2.34 for j = 1, k = 4 → 64;
- 0.55 → 2.11 for j = 2, k = 6 → 64;
- 0.52 at (16, 6) (critic).
The draft's "c(k, j) → 1, ≈ 0.8–0.9" was wrong.

**Prior art.** The non-adaptive case is classical: k ≥ 3 random entries per row need Θ(n) rows (Calkin 1997;
Balakin–Kolchin–Khokhlov 1991; Pittel–Sorkin 2016). Croot–Granville–Pemantle–Tetali (Annals 2012, §8.2) show that
unlimited large primes cut the stopping time only by a constant factor. Prop TAIL's contribution is the adaptive,
entry-wise-conditional version with j forced entries. No prior statement of that version was found by us.

**Consequence (narrowed after review).**
- Under LQ and H_tail⁺, bounded forcing cannot reduce the number of relations needed for one dependency below
  2^{−8}/(kp) = L[β]^{1−o(1)}, the factor-base scale. So the collection term stays MNFS's.
- **Not established** (codex #9):
  - that every adaptive arithmetic strategy obeys H_tail⁺ (the joint law is an assumption, part of the obstruction);
  - forcing by other means (polynomial or field choice in MNFS);
  - any elimination-cost lower bound. "Peeling stalls at the 2-core" is heuristic for adaptive matrices.
- So open item (f) is reduced, conditionally, to: H_sel-type selection (forcing ≫ 1 band primes means index
  > U), sparse LA on L[β]^{1−o(1)}-dimensional matrices (heurnfs lever (b): at best 1.8595 at ω ≲ 1.67), and the
  randomness hypothesis H_tail⁺ itself.

---

## 4. PTE O2: pencils of two presentations (Task 4)

### 4.1 The naive three-term rigidity is false; a Lemma R-shaped version cannot bound O2 anyway

**Counterexamples.**
- (a) X(X+3)² − 2(X+1)³ + (X+2)(X−1)² = 0 [PROVED, by expansion; found by E1's pilot].
  - c = (1, −2, 1), G = 0, pairwise coprime, M = 3.
  - E1 finds 92 lines like it at (M, B) = (3, 8).
- (b) Codex: (X+8)(X+1)(X−3)(X−4) + (X+5)(X+4)X(X−7) − 2(X+7)(X+2)(X−1)(X−6) = −72 [PROVED, by expansion].
  - M = 4, pairwise coprime, not an exact pencil.
  - Translations X ↦ X − s give an infinite family with the same G.
  - E1 finds 7009 projected (constant-G) lines with e = M and a |c| ≤ 2 triple at (4, 6), 332 at (5, 5), 53 at
    (6, 4), 8 at (7, 4) and 2 at (8, 3), and none at M = 9, 10 (B = 3).
- All of these are outside a Lemma R-shaped regime (M ≥ 2⌊log₂|G|⌋ + 4). E1 finds **no** in-regime pairwise-coprime
  small-c triple (§5.1, after a disclosed rescoring).

**Why a three-term rigidity lemma is the wrong tool for O2 [PROVED, elementary].**
- O2 members P_ν = F₀ + νF₁ + r_ν give, for three members, the identity Σc_iP_i = Σc_ir_i with c ∝ ν-differences.
- So a rigidity lemma with |c_i| ≤ C can only forbid triples whose ν-ratios (ν₁ − ν₀):(ν₂ − ν₀) have height ≤ C.
- For any fixed C, ν-sets of size ν_max^{1−o(1)} avoid every such triple. Each such equation says one element is a
  weighted average of the other two, so one Behrend sphere with a large no-carry digit base avoids all of them at
  once (PROVED from KNOWN: Ruzsa 1993 Thm 2.3 covers each single equation of type (2, 1); the simultaneous
  version is our straightforward adaptation). The construction is Behrend-type,
  applied to the finitely many invariant equations Σc_iν_i = 0, Σc_i = 0 [KNOWN in substance: Behrend 1946; Ruzsa,
  Acta Arith. 65 (1993)].
- So TR is consistent with ν-sets of size ν_max^{1−o(1)}, and TR alone cannot bound O2. The right statement is a
  *count*. (This shows a limitation of triple exclusions; it does not construct split O2 members.)

**Lemma TL (three-term localisation) [PROVED].**
- *Setting.* P₀, P₁, P₂ monic of degree M, split, roots in [−B, B]; c_i ∈ Z \ {0}; G = Σc_iP_i ≠ 0 with
  deg G < K and ‖G‖ ≤ h; g = gcd, U_i = P_i/g, r = M − deg g; Δ = 2^{K−1}√K h B^{K−1}, L = ⌊log₂Δ⌋.
- (a) If a is a root of exactly two of the U_i, the third has ≤ L roots outside {a − 1, a + 1}.
- (b) Two such points at distance ≥ 3 for the same third polynomial force r ≤ 2L, a common factor of degree
  ≥ M − 2L.
- (c) If G = 0, the U_i are pairwise coprime.
- *Proof.*
  - G̃ = G/g ∈ Z[X], with M(G̃) ≤ √K h, so |G̃(a)| ≤ Δ (Lemma R, step 3).
  - At a: c_lU_l(a) = G̃(a) ≠ 0, so 1 ≤ ∏|a − b| ≤ Δ. Every factor is a nonzero integer, and each factor ≥ 2
    doubles the product.
  - (b) is disjointness of the two exceptional sets. (c): a common root of two U's is a root of the third.
- *Consequence.* Every three-term configuration is
  - (A) a common factor of degree ≥ M − 2L;
  - (B) localised: some U_l is supported, up to L roots, on two consecutive-but-one integers. These are CUN-like
    rows on few pieces;
  - or (C) pairwise coprime: three members of one pencil (if G = 0) or of one linear system F₀ + ZF₁ + Z[X]_{<K}.
- In E1, class A dominates at (3, 8) and (4, 6), where it is vacuous because r ≤ M ≤ 2L. Class C dominates from
  M = 5 on (§5.1).

### 4.2 Lemma PEN (split members of a pencil) [PROVED from KNOWN]

**Statement.**
- *Setting.* F₀ ∈ Z[X] monic of degree M, 0 ≠ F₁ ∈ Z[X] with deg F₁ < M. T = gcd(F₀, F₁) has degree M − e; write
  G_i = F_i/T.
- (A) For t ≠ t′ the moving parts G₀ + tG₁ and G₀ + t′G₁ have no common root. Every root a satisfies G₁(a) ≠ 0,
  and determines t = −G₀(a)/G₁(a).
- (B) At most 2e − 2 members have a repeated moving root.
- (C) If e ≥ 3, the number of t ∈ Z for which F₀ + tF₁ splits over Z with all roots in [−B, B] is at most
  (16/3)(c·B^{1/2}(log B)^κ + 1) + min(2e − 2, 2B + 1). This is uniform in the heights of F₀, F₁; the additive
  term depends on e. Here c, κ are the
  absolute constants of Binyamini–Cluckers–Novikov.
- (D) e ≤ 2 is sharp: (X − a)·T and T·(X − a)(X − (s − a)) give ≥ B + 1 members.

**Proof.**
- (A) A common root α gives (t − t′)G₁(α) = 0, so G₀(α) = 0, contradicting coprimality.
- (B) Repeated roots occur only at critical values of −G₀/G₁. Equivalently, disc_X(G₀ + tG₁) is a nonzero
  polynomial in t of degree ≤ 2e − 2.
- (C) The ordered root pairs of a squarefree member lie on C: (G₀(x)G₁(y) − G₀(y)G₁(x))/(x − y) = 0, a curve of
  degree ≤ 2e − 2. Distinct members give disjoint pairs, by (A).
  - Lines on C:
    - vertical or horizontal lines would need a common root, or G₀ ∝ G₁;
    - the diagonal is excluded by the Wronskian G₀′G₁ − G₀G₁′ ≢ 0;
    - y = λx + ν forces G₀/G₁ to be invariant under x ↦ λx + ν. The stabiliser is finite, and over Q contains at
      most one involution x ↦ ν − x (two would compose to a translation).
  - So at most one rational line, carrying ≤ e pairs per member. Each squarefree member therefore has ≥ e(e − 2)
    pairs on the other components.
  - Rational components of degree δ_i ≥ 2 carry ≤ cδ_i²B^{1/δ_i}(log B)^κ integral points (BCN, arXiv 2302.04209,
    Thm 2). The others carry ≤ δ_i² (Bézout with a conjugate). Σδ_i ≤ 2e − 2.
  - Divide by e(e − 2) and use 4(e − 1)²/(e(e − 2)) ≤ 16/3. ∎
- Codex's blind round proves the same theorem independently (M ≥ 4, coprime pencils), with a sharper treatment of
  repeated roots: Σ_t(M − r_t) ≤ 2M − 2 via the Wronskian. This generalises ptebarrier's k = 1 fibre bound
  (F₁ = 1) to arbitrary pencils.

### 4.3 Prop O2-count (height-free) [PROVED from KNOWN]

**Statement.**
- *Setting.* F₀ monic of degree M, deg F₁ < M, k ≥ 0, and 𝓢 = {P ∈ F₀ + ZF₁ + Z[X]_{<k} split over Z, roots in
  [−B, B]}.
- For P ∈ 𝓢 let W be the product of its k largest roots, and U = P/W.
- Then U lies in the pencil 𝒫_W = {quot(F₀, W) + ν·quot(F₁, W)}, and P ↦ (W, U) is injective. Hence

  |𝓢| ≤ Σ_W N(𝒫_W) ≤ C(2B+k, k) · max_W N(𝒫_W),

  with:
  - N(𝒫_W) ≤ 1 if quot(F₁, W) = 0;
  - N(𝒫_W) ≤ 2B + 1 always, by PEN (A);
  - N(𝒫_W) ≤ PEN (C) if the moving degree e_W ≥ 3.

**Proof.** The quotient by a monic W of degree k depends only on the coefficients of degree ≥ k, and is linear.
The remainder r (deg < k) has quotient 0, and quot(P, W) = U because W | P. ∎

**Reading.**
- One pair of presentations with remainders of degree < k yields at most C(2B+k, k)·O(B^{1/2}(log B)^κ + M)
  non-degenerate members.
- This holds whatever the heights of F₀, F₁ and whatever the size of the remainders.
- Codex's (22) has a factor (2h+1)^k in the remainder height h instead. At k = 1, h can be e^{O(ℓ)} ≫ B, so
  removing it matters.
- At k = 1 the *nondegenerate* part (W-pencils with e_W ≥ 3) has ≤ (2B + 1)·O(B^{1/2}(log B)^κ + min(2e, 2B))
  members, which is B^{3/2+o(1)} when M = B^{o(1)}.
- Degenerate W-pencils (e_W ≤ 2) contribute ≤ 2B + 1 each, so up to O(B²) in total, and that is attained. Codex's
  example: M = 3, k = 1, F₀ = X³, F₁ = X; every root triple with sum 0 is a member, Θ(B²) of them.
- Hypotheses: 0 ≤ k ≤ M and W = ∏(X − a) over the k largest roots.
- **This does not refute O2.** L[β] = B^{1+o(1)} members are not excluded for k ≥ 1. It is *bounded amplification*,
  like ptesplit's mechanism B: ≤ B^{1/2+o(1)} members per top-root class W. Whether the pencils 𝒫_W can be made
  simultaneously rich remains open; E1 (§5.1) finds them poor.

---

## 5. Experiments (pre-registered, `notes/pilist_prereg.md`, commit c78c6bc)

Resources: E1 210 s, peak RSS 83 MB; E2 161 s, 170 MB; E3 (see §5.3). Memory pressure was normal (≥ 60% free)
before each run. The runs were sequential.

**Disclosed deviation (E1).**
- The first E1 launch crashed in the kproj = 1 branch, before writing any output. Members sharing a projection
  (same e₁..e_{M−1}, different constant) gave repeated line positions t and a zero gcd.
- Such triples are two-term k = 1 fibre relations, not three-term ones.
- Fix, applied before the rerun: in kproj = 1, equal projections are collapsed into classes, lines are found among
  distinct projections, and a triple must use three distinct classes. `small_c_triples` skips repeated t.
- kproj = 0 output (printed for (3, 8) before the crash) was identical in the rerun.

### 5.1 E1 — pencils and three-term identities

| (M, B) | split polys | lines kproj=0 / 1 | P1 violations | max members, e ≥ 3 (kproj=0 / 1) | e = M small-c lines (kproj=0 / 1) | TL classes of small-c triples, kproj=1 (A / B / C) |
|---|---|---|---|---|---|---|
| (3, 8) | 969 | 10218 / 24737 | 0 | 3 / 41 | 92 / 15462 | 53856 / 260 / 12188 |
| (4, 6) | 1820 | 16458 / 33291 | 0 | 3 / 28 | 0 / 7009 | 8526 / 143 / 20267 |
| (5, 5) | 3003 | 24013 / 24993 | 0 | 3 / 7 | 0 / 332 | 348 / 0 / 30173 |
| (6, 4) | 3003 | 14982 / 15142 | 0 | 3 / 5 | 0 / 53 | 51 / 0 / 23067 |
| (7, 4) | 6435 | 39534 / 39546 | 0 | 3 / 3 | 0 / 8 | 8 / 0 / 54879 |
| (8, 3) | 3003 | 12342 / 12345 | 0 | 3 / 3 | 0 / 2 | 2 / 0 / 16830 |
| (9, 3) | 5005 | 23067 / 23067 | 0 | 3 / 3 | 0 / 0 | 0 / 0 / 29799 |
| (10, 3) | 8008 | 40612 / 40612 | 0 | 3 / 3 | 0 / 0 | 0 / 0 / 50193 |

(Class-C counts in kproj = 1 include the G = 0 triples of kproj = 0: the pencils.)

**Evaluation.**
- **P1 held**: 0 violations of Claim A and of the 2e − 2 bound.
- **P2 held**: e ≤ 2 lines with ≥ B + 1 members exist in every cell (the reported maximum is capped at 9: 9, 9, 9,
  9, 9, 7, 7, 7 against B + 1 = 9, 7, 6, 5, 5, 4, 4, 4).
- **P3 FAILED at M = 4.**
  - It held at M = 3 (92 lines).
  - At M = 4 there is no *exact* coprime pencil (kproj = 0) with e = M, ≥ 3 members and a |c| ≤ 2 triple.
  - The weaker kproj = 1 form (constant G, codex's example (b)) holds at M = 4 with 7009 lines.
  - So the naive *exact* (G = 0) three-term statement fails at M = 3 only in our boxes. Its constant-G form fails
    up to M = 8.
- **P4 FAILED as scored by the registered code; post-hoc rescoring removes both flags.**
  - The code flagged 2 in-regime class-C triples at (6, 4), kproj = 1. Both have G = 0, e.g.
    {−3,−1,−1,2,2,4}, {−2,−1,−1,2,2,3}, {−1,−1,0,1,2,2} with c = (−1, 2, −1).
  - The registered in-regime test for G = 0 used the moving degree of the whole line (≥ 6). These *triples* share
    the fixed part (X+1)²(X−2)² of degree 4, so their own moving degree is 2. They belong to the sharp e ≤ 2
    family (Lemma PEN (D)).
  - With the triple's own moving degree (post hoc) there is no in-regime class-C small-c triple in any cell.
  - Every constant-G small-c triple with M ≥ 6 is out of regime or of class A.
- **P5 held**: max members on an e ≥ 3 exact pencil is 3 in *every* cell. In these boxes, nondegenerate exact
  pencils never carry more than 3 split members, and from M = 4 on none with e = M has a small-c triple. Lower e
  occurs: multiplying the cubic identity by (X − 2)^{M−3} gives an e = 3 pencil with c = (1, −2, 1) in every box
  (codex).

**Reading.**
- The naive three-term rigidity is false at small M, in both its exact and its constant-G forms.
- No in-regime counterexample exists up to M = 10 (roots in [−3, 3] at M ≥ 8).
- Exact nondegenerate pencils are very poor: ≤ 3 members, far below the PEN bound. This is consistent with an
  unproven conjecture that they carry O(1) members.

### 5.2 E2 — Johnson-radius decoders on the c60 rational side

Box smooth rate per coprime point: 0.0204 (X = 2^10), 0.0200 (2^13), 0.0253 (2^16).

| X | z | log P | log S | pred out/call | obs out/call | truth / pred (all points) | smooth outputs | Γ̂ | JR bound (prereg) | decoder misses |
|---|---|---|---|---|---|---|---|---|---|---|
| 2^10 | 13 | 10.3 | 8.87 | 0.273 | 0.230 | 46 / 54.6 | 2 | 0.49 | 8.9 | 0 |
| 2^10 | 20 | 16.1 | 11.08 | 0.267 | 0.255 | 51 / 53.5 | 1 | 0.25 | 11.3 | 0 |
| 2^10 | 30 | 22.6 | 13.12 | 0.161 | 0.170 | 34 / 32.3 | 4 | 0.98 | 9.8 | 0 |
| 2^10 | 60 | 45.0 | 18.53 | 0.0261 | 0.020 | 4 / 5.2 | 1 | 0.25 | 4.1 | 0 |
| 2^10 | 100 | 79.8 | 24.66 | 0.0026 | 0.005 | 1 / 0.53 | 0 | 0 | 1.07 | 0 |
| 2^10 | 160 | 137.8 | 32.41 | 1.3e−4 | 0 | 0 / 0.03 | 0 | 0 | 0.17 | 0 |
| 2^13 | 13 | 10.3 | 10.00 | 0.546 | 0.550 | 110 / 109.1 | 3 | 0.75 | 7.9 | 0 |
| 2^13 | 20 | 16.1 | 12.49 | 0.324 | **0.770** | **154 / 64.9** | 2 | 0.50 | 15.5 | 0 |
| 2^13 | 30 | 22.6 | 14.81 | 0.169 | 0.220 | 44 / 33.8 | 2 | 0.50 | 14.3 | 0 |
| 2^13 | 60 | 45.0 | 20.91 | 0.0229 | 0.015 | 3 / 4.6 | 2 | 0.50 | 5.3 | 0 |
| 2^13 | 100 | 79.8 | 27.82 | 0.0016 | 0 | 0 / 0.31 | 0 | 0 | 1.07 | 0 |
| 2^13 | 160 | 137.8 | 36.56 | 5e−5 | 0 | 0 / 0.01 | 0 | 0 | 0.12 | 0 |
| 2^16 | 20 | 16.1 | 13.77 | 0.270 | 0.281 | 18 / 17.3 | 0 | 0 | 19.2 | 0 |
| 2^16 | 30 | 22.6 | 16.32 | 0.212 | 0.172 | 11 / 13.6 | 0 | 0 | 20.8 | 0 |
| 2^16 | 60 | 45.0 | 23.04 | 0.0221 | 0 | 0 / 1.4 | 0 | 0 | 7.2 | 0 |
| 2^16 | 100 | 79.8 | 30.66 | 0.0011 | 0 | 0 / 0.07 | 0 | 0 | 1.15 | 0 |
| 2^16 | 160 | 137.8 | 40.29 | 2e−5 | 0 | 0 / 0.001 | 0 | 0 | 0.09 | 0 |

(2^16, 13) was skipped as registered (log S > log P). The decoder's certified threshold was within
[−4.7, +1.14] nats of log S. At z = 160 the ladder stopped at m = 16, 1.14 nats above the threshold. The 0-miss
check used the certified threshold.

**Evaluation.**
- **P6 held**: 0 misses in 272 decoder calls.
- **P7 FAILED in 1 of 8 eligible cells**: (2^13, 20), 154 against 64.9, ratio 2.37. The other seven ratios are
  0.84, 0.95, 1.05, 1.01, 1.30, 1.04 and 0.81.
- **P8 held**: max 0.77 outputs per call.
- **P9 held**: Γ̂ ≤ 0.98 in every cell, below each JR bound and below 25. The JR column above is the
  pre-registered table. The runner recomputed λ from the actual largest b and got slightly smaller values (e.g.
  8.35 vs 8.9); the table values were used for scoring.

**Post-hoc diagnosis of the P7 failure** (`src/pilist/posthoc.py`; `experiments/raw/pilist_E2_posthoc.json`,
`pilist_E2_coprime_posthoc.json`).
- The registered truth count included *non-primitive* points. In that cell 68 + 18 + 17 + 12 + 12 of the 154
  outputs lie on five divisor lattices whose shortest vectors are multiples of (30, 1): (−60, −2), (−300, −10),
  (−150, −5), (210, 7), (−420, −14).
- G(30, 1) = 30Y₁ + Y₀ has a large P-part, and the multiples k·(30, 1) inherit it. So the excess comes from
  non-primitive points. This diagnosis removes the specific discrepancy; it does not prove the conditional model.
- Restricted to gcd(a, b) = 1, with the primitive (projective) density Pr[q | G] = 1/(q + 1), the cells
  with pred_coprime ≥ 10 give observed/predicted 0.84, 0.64, 1.25, 1.12, 0.59, 1.25 and 0.40.
  - Total: 85 against 96.2, 0.88. The z-cells share points, so this is not a sum of independent counts.
  - Eligibility was redefined post hoc. It drops (2^16, 30), and the last retained primitive prediction is only
    4.98.
  - The density is exact on average over uniform primitive pairs. Per row it is 1/q for q ∤ b and 0 for q | b,
    which consecutive rows average correctly.
- The 1/(q + 1) density only lowers the probabilities, so Prop JR's hypothesis Pr[v_q ≥ j] ≤ q^{−j} is
  conservative for primitive points.
- The registered failure stands.

**Reading.**
- At Boneh's threshold a call has ≤ 0.77 certified outputs (all points) and ≤ 0.115 primitive ones (23/200 at
  (2^13, 13)), and ≪ 1 once P has more than about 10 primes.
- Pooled: 17 smooth certified outputs against 56.6 if every call were worth one box point, so Γ̂ ≈ 0.30. The
  per-cell Γ̂ rest on 0–4 events.
- The yield uses the threshold truth set (harness). The operational HG decoder was checked separately (P6).
- This is the toy-scale image of Prop JR: the call is a batched set of lattices of index ≥ V, each holding
  ≤ V/d ≤ 1 expected points.
- Toy scale only (u ≈ 3.5). The asymptotic claim rests on the proof.

### 5.3 E3 — adaptive closure of Θ(u) tails

Resources: 3692 s; peak RSS 72 MB. 20 seeds per cell; every run reached a dependency.

| (k, j) | strategy | median R*/n at n = 1000 / 2000 / 4000 / 8000 | runs with R* ≤ R_lb (of 80) |
|---|---|---|---|
| (3, 0) | nonadaptive | 0.917 / 0.919 / 0.919 / 0.920 | 1 |
| (5, 0) | nonadaptive | 0.992 / 0.992 / 0.993 / 0.993 | 0 |
| (5, 1) | nonadaptive | 0.992 / 0.992 / 0.993 / 0.993 | 0 |
| (5, 1) | greedy | 0.994 / 0.993 / 0.993 / 0.993 | 0 |
| (5, 1) | oracle (ILLEGAL) | 0.868 / 0.867 / 0.869 / 0.869 | 0 |
| (8, 0) | nonadaptive | 0.999 / 1.000 / 1.000 / 0.999 | 0 |
| (8, 2) | nonadaptive | 0.999 / 1.000 / 1.000 / 0.999 | 0 |
| (8, 2) | greedy | 1.000 / 1.000 / 0.999 / 1.000 | 0 |
| (8, 2) | oracle (ILLEGAL) | 0.892 / 0.896 / 0.897 / 0.898 | 0 |
| (12, 2) | nonadaptive | 0.999 / 1.000 / 1.000 / 1.000 | 0 |
| (12, 2) | greedy | 1.001 / 1.000 / 1.000 / 1.000 | 0 |
| (12, 2) | oracle (ILLEGAL) | 0.947 / 0.952 / 0.952 / 0.953 | 0 |
| (12, 4) | nonadaptive | 0.999 / 1.000 / 1.000 / 1.000 | 0 |
| (12, 4) | greedy | 1.000 / 1.000 / 1.000 / 1.000 | 0 |
| (12, 4) | oracle (ILLEGAL) | 0.827 / 0.834 / 0.835 / 0.836 | 0 |
| (16, 6) | nonadaptive | 1.000 / 1.000 / 1.000 / 1.000 | 0 |
| (16, 6) | greedy | 0.999 / 1.000 / 1.000 / 1.000 | 0 |
| (16, 6) | oracle (ILLEGAL) | 0.787 / 0.790 / 0.790 / 0.791 | 0 |

(In the j > 0 rows, `nonadaptive` is the j = 0 baseline at the same k, scored against the (k, j) threshold R_lb.)

**Evaluation.**
- **P10 held**: one legal run of 1040 had R* ≤ R_lb, at (2000, 3, 0) with R* = 1149 against 1425. At most 1 per
  cell, against the allowed 3.
- **P11 held**: legal medians R*/n ≥ 0.917. They are ≈ 0.92 for k = 3, the k = 3 XORSAT threshold region, and
  ≥ 0.99 for k ≥ 5.
- **P12 (exploratory)**: the ILLEGAL oracle, with forced entries chosen after seeing the random part, lowers the
  median only to 0.79–0.95. Legal greedy closing changes nothing measurable (0.993–1.000).

**Reading.**
- With bounded forcing, dependencies appear at R ≈ n, the column count.
- The proved threshold (R ≥ 2^{−8}/(kp); the exact first moment gives ≈ (0.5–2.7)/(kp)) is conservative here.
  The proved 2^{−8}/(kp) = n/(256k) is off by a factor of about 256k ≈ 700–4100; the exact first moment by about
  1.3–31.
- The j > 0 `nonadaptive` rows rerun the j = 0 protocol with the same seeds, so 240 of the 1040 legal runs are
  duplicates (≈ 800 distinct).
- Toy scale with uniform columns. H_tail⁺ is an assumption and is not tested against real smooth values.

## 6. Where the L[1/3] limit bites (synthesis; an interpretation, not a proved identification)

1. **Π_β itself (congruence form).** As a target, even beating MNFS requires the per-output cost to fall from the
   box sieve's L[1 − β] to L[1/3, O(1)]. That is essentially full output sensitivity. Only NFS-shaped listers admit
   fixed-fraction gains (Table 1.1).
2. **Box access in general.** The bite is a selection hypothesis. H_sel^gain is exactly what LQ needs, and an
   H_sel^gain violation with fresh primitive outputs, i.e. a violation of H_sel by ≥ (M/U)·L[θ′r]^{−1+ε}, would be
   an improvement (Remark SEL). We do *not* claim the two are equivalent.
3. **Certificate decoders.**
   - Conditionally on random divisibility on oblivious regions, for primitive points and in NFS scaling, the
     Johnson radius is where the limit bites for the
     whole Coppersmith/Boneh/GSS class: certified ⇒ d ≥ V, and batching gains L^{o(1)} (Prop JR; E2 Γ̂ ≈ 0.3).
   - The range of the method is proved optimal only for the univariate full-modulus case (Chinburg et al.).
4. **Adaptive Θ(u) tails.** Conditionally on H_tail⁺, the relation count cannot shrink (Prop TAIL). What remains
   is sparse LA with ω < 2 (heurnfs (b)), H_sel-type forcing, and the randomness hypothesis itself.
5. **PTE O2.** Counts, not rigidity: Lemma PEN and Prop O2-count give bounded amplification. Measured exact
   pencils carry ≤ 3 members.

Items 2–4 share one *form*: certifying a known factor-base divisor larger than the region searched should cost
about one unit per region-sized lattice. Each result proves this form only inside its own model or method class.

## 7. References

**Lattice and certificate decoding.**
- Coppersmith, J. Cryptology 10(4) (1997).
- Howgrave-Graham, IMA Cryptography & Coding 1997, LNCS 1355.
- Boneh, "Finding smooth integers in short intervals using CRT decoding", STOC 2000; JCSS 64 (2002) 768–784.
  Numbering is the JCSS version's: Thms 2.1, 3.1, 3.2, 4.1, §4.
- Goldreich–Ron–Sudan, "Chinese remaindering with errors", STOC 1999; IEEE TIT 46(4) (2000).
- Guruswami–Sahai–Sudan, "Soft-decision decoding of Chinese remainder codes", FOCS 2000, Thm 4.
- Guruswami, *List decoding of error-correcting codes*, PhD thesis, MIT 2001, §12.3.2. The section number is per
  codex (smoothenum_2); the literature critic could not retrieve the thesis.
- Cohn–Heninger:
  - "Approximate common divisors via lattices", ANTS X, Open Book Series 1 (2013);
  - "Ideal forms of Coppersmith's theorem and Guruswami–Sudan list decoding", AMC 9 (2015), Thms 1.1, 1.3.
- Herrmann–May, ASIACRYPT 2008, LNCS 5350, 406–424, Thm 4 (Assumption 1; polynomial time needs n = O(log log N)).
- Jochemsz–May, ASIACRYPT 2006, LNCS 4284, 267–282 (Assumption 1).
- Chinburg–Hemenway–Heninger–Scherr, "Cryptographic applications of capacity theory: on the optimality of
  Coppersmith's method for univariate polynomials", ASIACRYPT 2016 (arXiv 1605.08065), Thms 2, 13.
- Guruswami–Rudra, IEEE TIT 54(1) (2008) 135–150.
- Kopparty, "List-decoding multiplicity codes", Theory of Computing 11 (2015), art. 5.
- Brakensiek–Chen–Putterman–Zhang–Zheng, "Algorithmic list decoding of Reed–Solomon codes up to capacity",
  arXiv 2609.08005 (2026, preprint).
- Jeronimo, arXiv 2609.05870 (2026, preprint).

**Hardness and limits of decoding.**
- Guruswami–Vardy, IEEE TIT 51(7) (2005).
- Cheng–Wan: FOCS 2004 / SICOMP 37 (2007), the discrete-log reduction; ICALP 2008 / IEEE TIT 56 (2010), bounded
  distance at (2/3 + ε)·distance.
- Gandikota–Ghazi–Grigorescu, "NP-hardness of Reed–Solomon decoding and the Prouhet–Tarry–Escott problem", FOCS
  2016, 760–769; SICOMP 47(4) (2018) 1547–1584 (verified by the literature critic's second pass). The proof goes
  through Moments Subset Sum, and GGG identify PTE solutions as the obstacle to widening the gap: a thematic link
  to §4.
- Ben-Sasson–Kopparty–Radhakrishnan, IEEE TIT 56(1) (2010) 113–120: over characteristic-2 fields with K = N^δ,
  superpolynomial lists already at agreement ≈ N^{√δ}, beyond Johnson.

**Enumeration lower bounds.**
- Aono–Nguyen–Seito–Shikata, "Lower bounds on lattice enumeration with extreme pruning", CRYPTO 2018, LNCS
  10992, 608–637. The bounds hold in the Gaussian-heuristic model.
- Hanrot–Stehlé, "Worst-case Hermite–Korkine–Zolotarev reduced lattice bases", arXiv 0801.3331 (2008). This is
  the tightness result; their CRYPTO 2007 paper is the upper bound.

**Rational points and separated-variable curves.**
- Binyamini–Cluckers–Novikov, IMRN 2024(11) 9256–9265 (arXiv 2302.04209), Thm 2, with c = c(2).
- Context (pointers from the literature critic, not checked by us):
  - Fried, Crelle 1973;
  - Bilu–Tichy, Acta Arith. 95 (2000), on f(x) = g(y);
  - Avanzi–Zannier, Compositio 2003;
  - S. D. Cohen, Proc. LMS 43 (1981): O(H^{1/2} log H) specialisations with a smaller Galois group, by parameter
    height and with polynomial-dependent constants.
  - Nearest B^{1/2}-type result: Cluckers–Dèbes–Hendel–Nguyen–Vermeulen, Forum Math. Sigma 13 (2025) e153
    (arXiv 2311.16871), Thm 1.10, counts reducible specialisations as ≪ 2^{118(d+1)^n}(log‖F‖ + 1)^{10}B^{1/2}.
    That is polylog in height but exponential in degree, and by parameter height.
  - Also Walkowiak (Acta Arith. 2005), Castillo–Dietmann (Acta Arith. 180, 2017) and Dèbes (Manuscripta Math.
    1996; integral-point method).
  - (Pointers from the literature critic; not read by us.)
  - Lemma PEN counts by root range, uniformly in height and (apart from the additive term) in degree. That is the
    difference from the results above; no prior statement of this form was found by us or the literature critic.

**Additive combinatorics.** Behrend, PNAS 32(12) (1946) 331–332; Ruzsa, Acta Arith. 65(3) (1993) 259–282,
Thm 2.3.

**Coppersmith and smooth values, further.** Bernstein, "Reducing lattice bases to find small-height values of
univariate polynomials", MSRI Publ. 44 (2008), generalising Boneh §4; Aono, ePrint 2012/108, on the optimality of
Coppersmith lattices. (Pointers; not read by us.)

**Smooth numbers and large primes.**
- de Bruijn, Indag. Math. 13 (1951); Hildebrand–Tenenbaum, JTNB 5 (1993);
- Croot–Granville–Pemantle–Tetali, Annals 175(3) (2012);
- Lenstra–Manasse, EUROCRYPT 1990 / Math. Comp. 1994;
- Leyland–Lenstra–Dodson–Muffett–Wagstaff, ANTS V (2002).

**Random sparse F₂ systems** (context and non-adaptive baseline for Prop TAIL).
- Calkin, "Dependent sets of constant weight binary vectors", CPC 6(3) (1997) 263–271: random weight-k rows are
  independent for m ≤ βn, β < β_k, with β_k ~ 1 − e^{−k}/ln 2;
- Balakin–Kolchin–Khokhlov, Diskret. Mat. 3(3) (1991);
- Darling–Penrose–Wade–Zabell, EJP 19 (2014), art. 82;
- Kolchin, *Random Graphs*;
- Dubois–Mandler, FOCS 2002;
- Pittel–Sorkin, CPC 25(2) (2016);
- Dietzfelbinger et al., ICALP 2010, LNCS 6198.
- Prop TAIL's adaptive, forced-entry version goes beyond these non-adaptive results.

**Internal.** notes/smoothenum_results.md (LQ, BC, T1–T8, §3.7); notes/nonexact_results.md §6;
notes/ptecons_results.md §7; notes/ptebarrier_results.md (Lemma R, O2); notes/ptesplit_results.md;
notes/heurnfs_results.md.

## 8. Reviews

### 8.1 Codex blind round (`codex/pilist_1.out`, 18 min)

Same verdict; probability about 5%.

**Adopted:**
- the pencil theorem, proved independently (M ≥ 4, coprime; sharper repeated-root count);
- the M = 4 constant-remainder counterexample and its translation family;
- the first-moment tail theorem in support-counting form;
- the hidden-residue and public-m remarks on H_sel;
- the output-capacity remarks;
- the H_D calibration at α = 1/4 (eqs. 8–10);
- the block-folding estimate (eq. 17);
- the Chinburg–Hemenway–Heninger–Scherr and Herrmann–May citations.

**Ours:**
- the canonical binary-form statement and its congruence/NFS split;
- Lemma CS and the Johnson-radius ⇒ d ≥ V mechanism;
- the regime proof of unweighted JR, now superseded;
- Remark SEL;
- Prop O2-count (height-free) and Lemma TL;
- the Behrend limitation;
- the explicit tail constant;
- all experiments.

### 8.2 Codex adversarial round (`codex/pilist_2.out`): MAJOR REVISION — all points applied

| # | Severity | Finding | Fix |
|---|---|---|---|
| 1 | FATAL (to JR as stated) | supersets permit Γ = V | JR bounds certified outputs only |
| 2 | MAJOR | divisor sum chained into the exact-gcd moment | moment bound and divisor majorant kept separate |
| 3 | MINOR | uniformity in the regime proof | superseded by the weighted proof |
| 4 | MAJOR | H_div does not cover arithmetic selection | oblivious regions; "proves H_sel" withdrawn |
| 5 | MINOR | the weighted proof; GSS tolerance term | **weighted JR proof adopted** (checked by us); GSS form corrected |
| 6 | FATAL (to EQ-SEL) | the equivalence fails: θ = 0 sieving also violates H_sel | Remark SEL replaces EQ-SEL |
| 7 | MAJOR | nonempty lattice without a primitive point (E = 2, M = 12) | primitive-point lemma restated |
| 8 | MAJOR | TAIL row lengths and entry-wise conditioning | H_tail⁺, exactly k |
| 9 | MAJOR | (f) not reduced unconditionally | conditional consequence |
| 10 | MINOR | O2-count degenerate O(B²) (X³ + νX example) | reading fixed |
| 11 | MAJOR | e ≤ 3 cap not enforced in `targets.py` | enforced; Table 1.1 corrected (0.8246 / 0.5732; θ* 0.0303 / 0.0221) |
| 12 | MINOR | conditional costs | labelled |
| 13 | MAJOR | E3 missing at review time | run (§5.3) |
| 14 | MINOR | e = M qualifier; harness yield | applied |
| 15 | MAJOR | the capacity-decoding literature claim was outdated (Brakensiek et al. 2026); Cohn–Heninger 1.3 misclassified | §2.3 rewritten |
| 16 | MINOR | Hanrot–Stehlé, Cheng–Wan and Kopparty citations; g_n ≤ μ² by AM–GM | applied |
| 17 | MAJOR | "one wall" is an interpretation | §6 retitled and narrowed |

### 8.3 Opus proof critic: MAJOR REVISION, no FATAL — all points applied

**MAJOR:**
- M1: H_div is circular for arbitrary regions; restricted to oblivious regions and primitive points.
- M2: Π-nonexact is outside JR's scope (counterexample in congruence scaling); removed.
- M3: EQ-SEL is overstated, and the primitive-point lemma was wrong. Replaced by Remark SEL with H_sel^gain.
- M4: Herrmann–May and Chinburg et al. scope narrowed; the homogeneous-form caveat added.
- M5: extra outputs. JR now bounds certified outputs only.
- M6: the draft was incomplete; now complete.

**MINOR:**
- m1: the TAIL constants were wrong; replaced by a proved c₀ = 2^{−8} and the measured table.
- m2: H_tail⁺ clarifications.
- m3: regime (iii) constants; superseded.
- m4: the gain type is exp(O(u + …)), not (u log u)^{v_box}.
- m5: degenerate W-pencils.
- m6: the Behrend wording.
- m7: θ* thresholds; numerical labels.
- m8: E2 reporting (0.115 primitive; eligibility disclosure; pooled Γ̂ ≈ 0.30).
- m9: the 1/(q+1) density is exact on average.
- m10: the class-A dominance claim.

**NITs:** applied, including g_n ≤ μ² upgraded to PROVED and V versus X up to a factor of 2.

The critic verified:
- the exponent algebra of all three draft regimes;
- Lemmas S≥V, CS, TL and PEN (A)–(D), and O2-count;
- both identities;
- the pre-registration timestamps;
- the E1/E2 scoring against the raw files.

### 8.4 Opus literature critic: OK WITH FIXES — all applied

- All cited works exist.
- The Brakensiek et al. preprint is confirmed, including its derivative-free claim. The draft's "lists of bounded
  size" was wrong (n^{O_ε(1)}, time q^{ε^{−O(1/θ)}}); fixed.
- §7 corrections: Hanrot–Stehlé, Cheng–Wan and Kopparty. The two passes disagreed on the GGG year; we use the
  verified FOCS 2016.
- Added:
  - Guruswami's thesis §12.3.2, Boneh Thm 3.2 and the Cohn–Heninger multi-input ACD (beyond Johnson only with
    correlated inputs);
  - Jeronimo 2026, Ruzsa Thm 2.3 (the simultaneous version is ours), Calkin 1997;
  - the large-prime cycle literature;
  - the Chinburg et al. quote that the divisor version is open.
- The JM extrapolation is labelled.
- Novelty: none claimed; "no prior statement found by us" is defensible for the JR impossibility and the
  height-free pencil count.
- The critic read §0 and §6 before their rewrite. The current versions contain none of the withdrawn claims.

### 8.5 Confirmation round (Opus proof critic, on the revised text): ACCEPT WITH FIXES — all applied

**Verified correct:**
- the weighted Prop JR proof, every step (every error term is O(u) or smaller);
- the Prop TAIL constant proof (bracket ≤ e·2^{−4−k}, total ≈ 0.04);
- Remark SEL, H_sel ⇒ H_sel^gain, and the corrected primitive-point lemma;
- the E3 table against the raw file.

**Minor fixes applied:**
- the class output is defined as Cert;
- the threshold for a useful H_sel violation is stated in §0 and §6;
- "saves no relations" → "saves at most an L^{o(1)} factor", with the full list of remaining escapes;
- the NFS-scaling and primitive scope added in §0 and §6;
- relations with fewer band primes are handled;
- the E3 conservative factor corrected (256k vs exact 1.3–31) and duplicate runs disclosed.

**Nits applied:** the absorbed error terms are named; π(M) is defined.

The literature critic's final pass (OK WITH FIXES) is also applied: the GSS-type extrapolation label, the GGG and
BKR wording, the Calkin gloss, the thesis section flagged as per codex, the non-adaptive prior art in §3.3, and the
CDHNV/Walkowiak context for PEN.

No further codex confirmation round was run on this revision.

## 9. Reproduction

```bash
PYTHONPATH=. .venv/bin/python -m pytest -q tests/test_pilist.py
PYTHONPATH=. .venv/bin/python -m src.pilist.targets              # Table 1.1 (~3 min)
PYTHONPATH=. .venv/bin/python -m src.pilist.run_experiment E1     # ~4 min, 83 MB
PYTHONPATH=. .venv/bin/python -m src.pilist.run_experiment E2     # ~3 min, 170 MB
PYTHONPATH=. .venv/bin/python -m src.pilist.run_experiment E3
PYTHONPATH=. .venv/bin/python -m src.pilist.posthoc               # POST HOC: E2 divisor diagnosis
# coprime_recount() and primitive_model_check() in src/pilist/posthoc.py (POST HOC)
```
