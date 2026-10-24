# uprime: pre-registration (E1, E2)

Agent `uprime`, phase 26, 2026-10-24. Written and committed BEFORE either experiment is run. Code is frozen at this
commit: `src/uprime/{hybrid_bv,toy,experiments,analyze}.py`. Only smoke tests on non-registered parameters (x = 2^22;
seed 999 / 4242) were run before this commit, to check the code.

Both experiments are **toy checks** of premises. They support nothing asymptotic. The hybrid restricted BV and the
U′ re-proofs are mathematical statements (notes/uprime_results.md); a toy cannot test an error term of size
exp(−cD^{3/8}/log²D).

## E1 — hybrid restricted BV at toy scale (`src/uprime/hybrid_bv.py`)

**Question.** Do moduli t = ℓ·s whose cofactor s is y₁-smooth (y₁ ≪ y₂) behave, for y₂-smooth values in residue
classes, like the full family of moduli ℓ·s of the same size? Lemma 8.5H (hybrid) says the error summed over the
y₁-smooth subfamily is small *relative to that subfamily's own density*. The danger it excludes is concentration of
bad moduli on smooth cofactors.

**Config.** x = 2^28, y₂ = 200 (values), y₁ ∈ {7, 13} (modulus cofactors), ℓ prime in (256, 512] (so ℓ > y₂),
Q ∈ {2^12, 2^13, 2^14} (Q ≤ √x), s ∈ (Q/2ℓ, Q/ℓ], t = ℓs ∈ (Q/2, Q]. All y₂-smooth n ≤ x are enumerated exactly.
Families per Q: F₁(y₁) = {s y₁-smooth}, F_all = all s (here s ≤ 64 ≤ y₂, so F_all is also the y₂-smooth-cofactor
family that Lemma 8.5′ at y = y₂ would control). Six cells (Q, y₁).

Per modulus: μ_t = Ψ_t/φ(t); maxE(t) = max over reduced r of |count(r) − μ_t|; dev = maxE/μ_t; bad = [dev > 1/2]
(paper2 def:good); V(t) = Σ_r (count − μ)²/(φ(t)μ); z(t) = maxE/√(2μ log φ(t)).
H(F) := Σ_{t∈F} maxE(t) / (Ψ(x, y₂)·#F/Q) (the hybrid sum normalised by the family's own density).

**Predictions (each must hold in all six cells).**
- P1a: no bad modulus in F₁.
- P1b: H(F₁)/H(F_all) ∈ [0.8, 1.25].
- P1c: median V over F₁ / median V over F_all ∈ [0.8, 1.25].
- P1d: median z over F₁ ∈ [0.3, 1.5].

Command: `.venv/bin/python -m src.uprime.hybrid_bv experiments/raw/uprime_e1.jsonl`, then
`.venv/bin/python -m src.uprime.analyze E1 experiments/raw/uprime_e1.jsonl`.

## E2 — mixed-smoothness acceptance and the conditional child mean (`src/uprime/experiments.py`)

**Question.** (a) In a toy design-(a‴) field, are rational B-smoothness and algebraic B₂-smoothness of level-1 box
points independent to first order (the product shape of FM-U′/LD*-U′)? (b) How much does the squarefree /
discriminant condition remove? (c) Is the mean number of good points in a child's box, *conditional on the parent being
accepted and q being its child*, close to the unconditional (product) main term measured by random labels of fresh
fields (Hyp CM, the conditional-mean bullet of U′)? (d) Are irregular child lattices rare? (e) Is the child count
Poisson-like at 0?

**Config.** 48-bit balanced semiprimes, seeds 1–6, 2 fields each (`toy.sample_field`, d = 3, M_bits = 20, H_c = 4),
300 special-q primes ℓ ∈ (2^17, 2^18] per field, level-1 reduced-basis box radius S₁ = 16, B = 2^11, B₂ = 2^16.
Acc: G/ℓ B-smooth, F B₂-smooth, the (B, B₂]-part of F squarefree and coprime to b·f_d·disc(f)·N. For the FIRST
accepted point per box (box order), every algebraic child (q, r): square child box of radius V₂ = ⌊√q·S₂⌋, S₂ = 24,
Good(P′): q ∥ F(P′), G(P′) and F(P′)/q B-smooth. Null: for each child, a fresh field from a per-seed pool of 4 and a
uniformly random simple label (q′, r′) with q′ ∈ [q/1.05, 1.05q], same box rule.

**Predictions.**
- P2a: pooled #both / Σ_boxes(#rat·#alg2/#points) ∈ [0.8, 1.3].
- P2b: pooled (removed_sq + removed_disc)/#both ≤ 0.02.
- P2c: R_CM := mean conditional child count / mean null child count ∈ [0.8, 1.5], and its 95% bootstrap CI (parent
  clusters) ⊂ [0.6, 2.0].
- P2d: fraction of children with an irregular lattice (4S₂²λ₁² < q) ≤ 3/S₂² ≈ 0.0052.
- P2e: Pr[conditional child count = 0] / e^{−mean} ∈ [0.8, 1.3].

Command: `.venv/bin/python -m src.uprime.experiments experiments/raw/uprime_e2.jsonl 6 1`, then
`.venv/bin/python -m src.uprime.analyze E2 experiments/raw/uprime_e2.jsonl`.

**Reading fixed in advance.** P2c is the only toy probe of Hyp CM. Holding, it says the parent's acceptance does not
deplete its children's boxes at toy scale; failing low (R_CM < 0.8) would be a warning that the conditional mean is
genuinely below the product main term, which the proof plan for CM would have to explain. Neither outcome bears on
the asymptotic statement. Deviations from this protocol will be disclosed in notes/uprime_results.md.
