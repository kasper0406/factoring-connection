# rank3_1 — blind round: rank r ≥ 3 linked modules between two nonlinear orders

You are a research number theorist (Diophantine approximation / geometry of numbers / decomposable forms).
You have read-only access to the repo at /Volumes/git/factoring-connection; scratch work goes in your cwd.
Web search is enabled; use it for literature. Honesty matters far more than a positive verdict: a sharp,
correct partial theorem or a precise statement of why the problem is out of reach is a valuable outcome.
Background notes you may read: notes/resbarrier_results.md (§2.3 Lemma C2, §2.5 Lemma KL, §2.6, §4.2),
notes/heurnfs_results.md (§1, §4). Do NOT read notes/rank3_* (another agent's working files).

## Setting (the "linked-order model", rank r)
- K, K′ number fields of degrees d, d′ ≥ 2 (both nonlinear), orders O ⊂ K, O′ ⊂ K′, ring maps
  φ: O → Z/N, φ′: O′ → Z/N (N squarefree, all prime factors > d + d′).
- Generators w_j = (u_j, u′_j) ∈ O × O′, j = 1..r, with φ(u_j) = φ′(u′_j) =: v_j ∈ Z/N ("linked"), and
  (v_1..v_r) primitive mod every p | N.
- Candidates: a ∈ Z^r in a box |a_j| ≤ E; θ(a) = Σ a_j u_j ∈ O, θ′(a) = Σ a_j u′_j ∈ O′.
- Smoothness integers: W(a) = N_K(θ(a)), W′(a) = N_{K′}(θ′(a)). These are decomposable forms of degrees
  d, d′ in r variables: W = ∏_σ ℓ_σ(a), ℓ_σ = σ∘θ, and W′ = ∏_τ ℓ′_τ(a). F := W·W′, degree D = d + d′.
- Cost model: candidate a is B-smooth with probability ρ(log|F(a)|/log B) (heuristic H); cost ≥ #examined
  + B². Write x = log N.

## What is known (in the repo)
- r = 2: every rank-2 linked module has N | Res(W, W′) hence N² | disc F; Bean (Compositio 1994, Thm 1)
  gives area{|F| ≤ 1} ≤ 16·|disc F|^{−1/(D(D−1))}, and Davenport/Thunder lattice-point counting turns this into
  #{box points with 0 < |F| ≤ X} ≤ 16·N^{−2/(D(D−1))}X^{2/D} + O(D·E). This yields the ledger
  log X ≥ x/(D−1) + (D/2)·log(#candidates) and α ≥ 1/3 (L[1/3, (32/9)^{1/3}] floor), rigorously in the model.
- Rank r with a *rational* side (d′ = 1): a Liouville/Minkowski argument on ker W′ gives a height bound
  (resbarrier Lemma KL), and r = 2 is optimal.
- Rank r ≥ 3 with both sides nonlinear: OPEN. For any rank-2 sublattice P with short basis, N | Res(W|_P, W′|_P),
  so a height bound d′h + d h′ ≥ 1 holds for the family, independent of r. This does not exclude an
  "aggregate" ledger log X ≈ x/d + (D/r)·C, which would give a formal L[1/4]. Only a Gaussian/Schmidt
  counting heuristic G_r (counts of degree-d algebraic points of bounded height + independence of the r − 1
  linking conditions) excludes it, and only for typical N.
- A second variant (from a different line, "split side"): the O-side norm is smooth *by construction*
  (θ(a) is a product of small linear pieces), so only W′ needs to be smooth; a family may then load
  height onto the O-side. G_r alone does not exclude this.

## Task
1. **Literature.** Find the exact theorems that count integer solutions of decomposable form inequalities
   |F(x)| ≤ X in r ≥ 3 variables, and that bound the volume V_F = vol{x ∈ R^r : |F(x)| ≤ 1}:
   Thunder, "Decomposable form inequalities", Ann. of Math. 153 (2001); Thunder's finite-type criterion;
   Evertse, Evertse–Győry (book "Discriminant equations" / "Unit equations"); Schmidt's norm-form theorems;
   any 2010–2026 work (e.g. Liu, Thunder–Liu, Stewart, Akhtari, Bennett–Dummigan–?...). State each theorem
   precisely: hypotheses, the invariant controlling V_F (is there a discriminant-type invariant for
   decomposable forms in r variables? which minors enter?), and the explicit dependence of counts on F.
2. **Ledger.** Decide whether an area/volume bound of the shape V_F ≤ c(r, D)·N^{−κ} with κ ≍ r/D² (or any
   bound strong enough to give α ≥ 1/3 for rank r ≥ 3 with both sides nonlinear, counting points in a box of
   side E ≤ B^{1−o(1)}) can be PROVED. The linking gives: for every u, v ∈ Z^r,
   N | ∏_{σ,τ} det(ℓ_σ, ℓ′_τ)(u ∧ v) (a decomposable form on Λ²). Can this be converted into a volume or
   counting bound? If not, find a counterexample (a linked rank-3 module, both sides nonlinear, with V_F
   large or many small values, at least in principle), or prove the strongest partial theorem: e.g. r = 3,
   fixed degrees, an exceptional set of N, a weaker exponent, or a bound on a different statistic
   (e.g. the number of candidates with BOTH |W| ≤ X and |W′| ≤ X′).
3. **Split side.** Same question when only W′ is charged (W is "free") — what is the right statement?
4. **Construction.** Is there any explicit rank-3 design (e.g. f, g of degree d with common root m mod N,
   candidates a + bα + cα² on both sides; or non-power-basis modules) whose (counting, smoothness)
   ledger beats rank 2 in the heurnfs model, even by a constant? Compute it.
Label every claim KNOWN (citation), PROVED (proof), HEURISTIC (assumption named) or CONJECTURE.
Be dimensionally careful (heights in units of log N, degrees d ≍ (x/log x)^{1/3}, etc.); you may run python.
Keep it under ~5000 words.
