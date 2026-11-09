# midk — pre-registration

Agent `midk`, phase 31 follow-up, 2026-11-09. Committed before any experiment is run or analysed.
Target (orchestrator): item #2, the **uncovered intermediate-k regime** — a regime where the
L[1/3] (or 2/5, 3/8) barrier rests on a counting/independence hypothesis at an INTERMEDIATE value of
a discrete parameter k, not on a proved impossibility. Identify it precisely; decide whether any
intermediate k can beat L[1/3].

Reused by import (never edited): `src/resbarrier/ledger.py` (G_r, KL, alpha_floor, gae/res floors),
`src/rankdegen/{fh,ledger}.py` (frame–Hölder, Φ program, LP_exc scope), `src/multidescent/model.py`
(α_j = (j+1)/(3j+2) ladder). Helpers I need from these are imported, or copied into `src/midk/`.

## Mapping (to be stated precisely in results, not pre-judged here)

Candidate intermediate-k barrier arguments whose only block at intermediate k is a
counting/independence hypothesis (as opposed to a proved area/discriminant/Coppersmith wall):
- **(A) rank-r linked modules (G_r regime)** — `resbarrier` §2.6/§4.2, `rankdegen` §5–6. Parameter
  k = rank r, intermediate 3 ≤ r ≤ d. Hypothesis: G_r (Schmidt point counts + independence of the
  r−1 linking congruences mod N); `rankdegen`'s LP_exc for large boxes; frame-less configs open.
- **(B) descent ladder** — `multidescent`/`pc3`/`c12rat`/`c12alg`. Parameter k = descent level j,
  intermediate 2 ≤ j < ∞. Hypothesis: pair-correlation C₁/C₂. Cost α_j = (j+1)/(3j+2).
- **(C) floor-sum ladder** — `floorhard`/`floorsum`/`topblock`. Parameter k = degree j / box side.
  Deterministic N^{1/(j+2)} metric (not L[1/3]); block is the Coppersmith known-bits barrier.

## Pre-registered predictions

**E1 (exact exponent ledger; `src/midk/ledger.py`, `experiments/raw/midk_e1.jsonl`).**
For each family, compute exactly, as a function of k, (i) the rigorous floor, (ii) the floor under
the counting hypothesis, (iii) the exponent if the hypothesis fails favourably. Predictions:
- **P1a.** Family A fallback (G_r fails; only the rigorous 2-plane resultant height floor x/d
  survives): the per-candidate log-norm obeys log X ≥ 2√(x·c), c = log(E√r) ≥ 1 (AM–GM on the
  height term x/d and the degree term d·c). Hence **α ≥ 1/4 rigorously**, and α = 1/4 is attained at
  d ≈ √x, r ≈ x^{1/4}(log x)^{3/4}. No value of (d, r) gives α < 1/4. (exact check of the AM–GM floor
  and of the attaining scaling; u = x^{1/4} at the optimum.)
- **P1b.** Family A under G_r (height floor (r−1)x/d): the optimum is at r = 2 and gives **α = 1/3**;
  α(r) is non-decreasing in r (√r growth). KL (rational side) and gr_min floors are non-decreasing
  in r. So **no intermediate r beats 1/3 under G_r**.
- **P1c.** Family B: α_j = (j+1)/(3j+2) is strictly decreasing in j with inf = 1/3 (j→∞) and
  α_j > 1/3 for every finite j; endpoint j=1 gives 2/5. So **no finite intermediate j reaches or
  beats 1/3** — favourable resolution of C₁/C₂ still gives α_j ∈ (1/3, 2/5).
- **P1d.** Family C is in the N^{1/(j+2)} deterministic metric; N^{1/(j+2)} ≫ L[1/3] for every fixed
  j, so it is not a route to beat L[1/3] (recorded, not a route).
- **Conclusion predicted:** the ONLY intermediate-k regime whose favourable failure yields α < 1/3
  is Family A (rank-r / G_r), with fallback exactly L[1/4].

**E2 (exact independence / construction test; `src/midk/independence.py`,
`experiments/raw/midk_e2.jsonl`).** Small N = pq, exact arithmetic (sympy/python-flint).
- **P2a (independence baseline).** For random linked rank-r modules over two degree-d number fields
  reduced mod a prime p, the fraction of modules whose r generators share a common nontrivial
  projective root (are "linked at rank r") is p^{-(r-1)}(1 + o(1)); over ≥ 10^4 trials the empirical
  count lies within a 3σ Poisson band of the independence prediction, for r = 3, 4 and several p.
- **P2b (structured dependence costs height or splits N).** For the structured constructions that
  §4.2 names as the only candidates to revive the family — (i) Galois/CM twin w′_j = σ(w_j),
  (ii) cyclotomic, (iii) "N is a norm" full-rank — each either (I) pays geometric family height
  (height ≥ N^{c/1} with c bounded below, measured exact), or (II) realises the dependence only
  through a Z/N-module endomorphism/automorphism that is non-scalar, in which case the registered
  routine `split_from_intertwiner` returns a nontrivial gcd with N (splits N). No construction gives
  dependence at height N^{1/d} without splitting N.
- **P2c (excess lies on sub-families).** When an engineered linked module shows more smooth-normable
  candidates than the rank-r volume bound predicts (excess over volume), the excess candidates lie on
  a proper sublattice of rank < r (so they are charged by a lower-rank ledger and add no new
  relation rank), confirming `rankdegen` §6.1–6.2 at intermediate r. Measured exactly on ≥ 8 toy
  linked modules; prediction: 100% of excess on sub-rank lattices, 0 genuine rank-r excess.

**Prop MID (candidate new statement; needs independent Opus critic before being reported as
proved).** A non-scalar Z/N-module endomorphism intertwining the two sides of a linked module
(the structure that would make the r−1 linking conditions dependent over the whole family) yields,
through its CRT eigen/idempotent decomposition, a nontrivial factor of N; equivalently, realising
such a dependence mod N is factoring-equivalent. Scope and hypotheses to be stated exactly in the
results note. Tested by E2 P2b(II).

## Status labels, arithmetic, reproducibility

Exact arithmetic only (sympy, python-flint, Fraction); floats for plots only. No leakage: the
independence/construction code receives N and public config; `split_from_intertwiner` is a
diagnostic (it may read p, q) and is named as such, never imported by any operational routine.
Raw data JSONL under `experiments/raw/`. Any claimed theorem passes a read-only Opus critic. Codex
GPT-6.1-sol: blind `codex/midk_1.md` (launched before analysis), adversarial `codex/midk_2.md`.
Commit locally by explicit path; do not push.
