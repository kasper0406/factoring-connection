# pc3 — pre-registration (experiments E1, E2)

Agent `pc3`, phase 26, 2026-10-24. Code: `src/pc3/toy.py`, `src/pc3/experiments.py` (frozen at the commit of this
file). Tests: `tests/test_pc3.py`. Raw output: `experiments/raw/pc3_e1.jsonl`, `experiments/raw/pc3_e2.jsonl`.

## What is tested

Only the RATIONAL block of multidescent's Hyp PC₃, at toy scale, with exact smoothness tests. These experiments test the
*shape* of the objects in the decoupling reduction (Prop DC, to be written in `notes/pc3_results.md`). They do not test
any asymptotic statement and support nothing rigorous. The algebraic block is not simulated.

Toy objects (`src/pc3/toy.py` docstring): g uniform on primitive vectors of [M/2, M]²; a label (q, r), q prime in
(500, 1000], r uniform mod q; reduced basis e₁, e₂ of L_𝔮 = {a ≡ rb mod q}; child box = primitive P′ = ie₁ + je₂ with
max(|i|, |j|) ≤ S modulo ±; child point good iff |G_g(P′)| is y-smooth; parents = primitive P = i₀e₁ + j₀e₂ with
K ≤ max(|i₀|, |j₀|) < 2K (a random sample of n_par, filtered so that |G| ≤ X); parent accepted iff |G_g(P)| = ℓ·s with ℓ
prime in (Y, 2Y] and s y-smooth. N(g) = number of good child points; Z(g) = number of accepted parents.

Frozen parameters: y = 150, Y = 2000, M = 2000, S = 20, Ks = (20, 80, 320), n_par = 1500, n_g = 300 000 (E1) and
60 000 (E2) independent g per label, X = 9·10⁷, 20 labels per experiment, seeds 20261024 (E1) and 20261025 (E2).

## Disclosure (before the main runs)

Smoke runs at seed 999 (1–2 labels, n_g = 20 000–40 000) were looked at while debugging. They showed that the first
local model (gcd(n₁, det) common to both values on the pair lattice, n uniform) misses a large negative covariance at
toy scale. Its source is the *primitivity of g*: for p ∤ det(P′, P″) the two values cannot both be divisible by p.
The model was therefore replaced by the **principal proxy**: 1_S(n) ↦ f(n) = τ_cls(small-prime valuation class of n)·
∏_{17≤p≤y} τ_p(v_p(n))·φ(|n|), with exact empirical local factors and the exact smooth density φ by size, evaluated on
the same sampled values. The bands below were set after seeing these smoke runs, so P1 and P2 are not blind.

## E1 — two-point structure of the child box (`python -m src.pc3.experiments E1 experiments/raw/pc3_e1.jsonl 20`)

For every pair of child points, C_emp = Pr(both good)/(Pr·Pr) − 1 over the g sample and C_mod is the same quantity for
the principal proxy. Pooled over the 20 labels (all pairs, equal weight per label); CIs from a label-cluster bootstrap
(2000 resamples).
- **P1a** mean residual ⟨C_emp − C_mod⟩ ∈ [−0.01, 0.01], with 95% CI ⊂ [−0.02, 0.02].
- **P1b** pooled slope cov(C_emp, C_mod)/var(C_mod) ∈ [0.8, 1.3].
- **P1c** excess dispersion var(C_emp) − noise − var(C_mod) ≤ 0.5·var(C_mod). Here noise = ⟨1/(n_g p′p″)⟩ is the
  binomial noise of C_emp.
- **P1d** box average ⟨C_emp⟩ ∈ [−0.01, 0.01]. Large per-pair covariances cancel over the box, so N is nearly Poisson.

## E2 — conditioning on the parent vs the unconditional child box (`… E2 experiments/raw/pc3_e2.jsonl 20`)

Per label: ϕ = Var N/m², ϕ_ex = ϕ − 1/m. Over parents with ≥ 30 accepting g:
- Φ = Σ_P Pr(E_P)·Var(N | E_P)/E[N | E_P]² / Σ_P Pr(E_P), and Φ_ex = Φ − ⟨1/E[N | E_P]⟩. This is the toy PC₃ quantity.
- shift = ⟨E[N | E_P]/m − 1⟩.
- σ² = Var Z/z², σ²_ex = σ² − 1/z.
- dc = E[Z(N − m)²]/(E Z·m²).

Pooled means over labels, for each K:
- **P2a** |⟨Φ_ex⟩ − ⟨ϕ_ex⟩| ≤ 0.01.
- **P2b** |⟨Φ_ex⟩(K = 320) − ⟨Φ_ex⟩(K = 20)| ≤ 0.01 (no growth with the parent's coefficient size).
- **P2c** ⟨σ²_ex⟩ ∈ [−0.02, 0.05] (the label's parent count is nearly Poisson).
- **P2d** |⟨shift⟩| ≤ 0.03.
- **P2e** |⟨dc⟩ − ⟨ϕ⟩| ≤ 0.01.

## Reading fixed in advance

If P2a/P2e hold, then at toy scale the weighted 3-point object of PC₃ equals the unconditional 2-point object up to the
stated tolerance, which is what Prop DC exploits. A failure would mean that toy-scale parent conditioning carries
covariance that the 2-point objects miss. Prop DC remains a theorem in either case; only its tightness is in question.
P1 tests whether non-principal contributions (characters, BDH/LSQ errors in the theory) are small next to the
principal local × archimedean covariance at toy scale. Toy u ≈ 2.9 and log y ≈ 5, so the local factors are large;
nothing here is asymptotic.
