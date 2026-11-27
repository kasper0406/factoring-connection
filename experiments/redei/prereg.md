# redei prereg (2026-11-27)

Agent `redei`. Goal: a poly(log N) native evaluator of a 2-power class-group bit of a discriminant attached to
N = pq that goes beyond genus, and a family that turns such bits into a factor.

## Hypotheses / predictions (registered before running)
- H1 (ledger, 4-rank). For D = -4kN, -8kN with k a product of public primes, r4 Cl(D) is a function of
  B0 = ((-1/p),(2/p),(p/q)) and x_r = (r/p), r | k. Prediction: over the family k in {1} ∪ {r prime <= B}
  the 4-ranks determine (r/p) for every r <= B up to the p<->q swap and a guess of B0 (O(1) bits), so the
  family = a Legendre oracle (./p) = a QR_N oracle (both directions, randomized). Bits: pi(B) + O(1).
  Kill for decoding: QR_N => factoring is an open classical problem; predict no decoder found.
- H2 (ledger, 8/16-rank). 8-ranks of the family add Redei-symbol bits [a,b,p] (QR in Z[sqrt a]/N) and
  p mod 16 in some cells (tower T3); predict 16-ranks do NOT add p mod 32 beyond what p mod 16 and
  lower symbols give (spin-type info instead). Measure: empirical conditional entropy H(p mod 32 | ...)
  on small N with PARI ground truth. Coppersmith needs p mod 2^j with 2^j > N^{1/4}; predict max j = O(1).
- H3 (native). Every route among {composition on public small-norm forms, square roots / 2-descent of
  public classes, Redei reciprocity / governing fields, Dirichlet/Hurwitz sums mod 2^k, Kronecker-Hurwitz /
  Eichler-Selberg relations mod 2^k} needs a square root mod N, a relation / group order, or >= N^c work.
  Each will be checked on small N against PARI and its native cost measured vs log N.
- Kill criterion (from brief): every route needs a sqrt mod N or a factorization-dependent object, or costs a
  power of N. Positive claim only with a factor-blind implementation given N alone on 35-bit+ inputs.

## Experiments
- E1: verify Redei-matrix 4-rank formula vs PARI quadclassunit on >= 300 (N, k) pairs, N < 10^7.
- E2: bit ledger for families F_B, B in {10, 30, 100}: count hidden-assignment preimages from 4-ranks.
- E3: 2-Sylow types (PARI) of Cl(-4kN), Cl(-8kN), k <= 30, on ~2000 small N: conditional entropy of
  p mod 2^j (j = 4,5,6), (p/q)_4, given public data + all 4-ranks, + all 8-ranks, + all 16-ranks.
- E4: native-cost probes: (a) public-form subgroup generation needs relations; (b) square root of the public
  class T via conic needs sqrt mod N; (c) h(-4N) mod 8 via Dirichlet sum (Gauss-lemma linearisation gives only
  mod 4); (d) Kronecker-Hurwitz sum mod 2^k: family size and per-term cost.
