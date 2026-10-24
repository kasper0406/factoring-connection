# Pre-registration: geometry (agent `geometry`, phase 26, 2026-10-24)

This note is committed **before** any code under `src/geometry/` is run on the instances below. The codex blind
prompt `codex/geometry_1.md` was launched before this note was written; its output had not been read when this
note was frozen.

## 0. Purpose

The phase-26 question is conceptual: which geometric features of objects built from N encode the factors, and can
they be computed without factoring? Most answers are proofs or literature (`notes/geometry_results.md`). The
experiments below are **correctness checks** of the formulas and reductions that the results note relies on. They
are not evidence of any speed-up. A reduction that factors N given an oracle value only says the oracle is as hard
as factoring (rule 8).

## 1. Instances

* `balanced_semiprime(n, seed)` from `src/common/instances.py`, n ∈ {32, 64}, seeds 0..199 (E2) and 0..1999 (E3).
* Small exact levels: all squarefree N ≤ 2000 (E1c), and all squarefree N with ψ(N) ≤ 400 (E1b).

## 2. Experiments

**E1. Modular-curve formulas (exact).**
* (a) Trace formula. `a_l(N)` := Tr(T_l | S_2(Γ_0(N))) for squarefree N, prime l ∤ N, computed from the
  Eichler–Selberg trace formula (Cohen's form, k = 2, trivial character):
  Tr T_l = σ(l) − (1/2) Σ_{t² < 4l} Σ_{f} h_w((t²−4l)/f²) μ(t,f,l) − (1/2) Σ_{dd'=l} min(d,d') Σ_{τ|N} φ(gcd(τ,N/τ)),
  with μ(t,f,l) = ψ(N)/ψ(N/N_f) · #{x mod N : x² − tx + l ≡ 0 (mod N·N_f)}, N_f = gcd(N, f), h_w the class number
  with weights 1/2, 1/3 at D = −4, −3. Checks: (i) N = 1 gives 0 for all primes l ≤ 200; (ii) N ∈ {11, 14, 15, 17,
  19, 21} against l + 1 − #E(F_l) for the Cremona curves 11a1, 14a1, 15a1, 17a1, 19a1, 21a1, all primes l ≤ 200,
  l ∤ N; (iii) N = 22 against 2·a_l(11a1); N = 26 against a_l(26a1) + a_l(26b1); N = 33 against a_l(33a1) + 2a_l(11a1).
* (b) Independent linear-algebra check. Build the weight-2 Manin-symbol space for Γ_0(N) (generators P¹(Z/N),
  relations x + xσ = 0, x + xτ + xτ² = 0) over F_P, P = 2^31 − 1, apply T_l by Merel's Heilbronn matrices
  {[[a,b],[c,d]] : ad − bc = l, a > b ≥ 0, d > c ≥ 0}, and compare (Tr(T_l) − (l+1)(2^{ω(N)} − 1))/2 with (a), for
  all squarefree N with ψ(N) ≤ 400 and primes l ≤ 23, l ∤ N. Also check dim = 2g + 2^{ω(N)} − 1.
* (c) Orbifold signature. For every squarefree N ≤ 2000 (and, as a bonus, every N ≤ 600): compute the permutations
  S: (c:d) ↦ (d:−c), ST and T on P¹(Z/N) explicitly; check |P¹(Z/N)| = ψ(N), #fix(S) = ν₂, #fix(ST) = ν₃,
  #cycles(T) = ν∞ against the closed formulas, and that Riemann–Hurwitz gives the genus formula.

**Prediction P1:** every check in E1 passes exactly (0 mismatches). If (a)(ii)/(iii) fails only for one curve, the
first suspect is a mis-remembered Cremona equation; we then check the discriminant's prime support and report the
discrepancy rather than editing the formula.

**E2. Reductions from geometric oracles to factors** (oracle values computed by a diagnostic harness from p, q; the
decoders receive only N and the oracle answers). n ∈ {32, 64}, 200 instances each. Oracles:
* G: genus g(X_0(N)); A: index ψ(N) (= 3/π · hyperbolic area of Γ_0(N)\H = |P¹(Z/N)|); P: Pillai
  Σ_{g=1}^{N} gcd(g, N) (total number of components of all star polygons {N/g}); L: σ(N) (number of index-N
  sublattices of Z²);
* K: necklace count oracle K(N, k, r) = Σ_{d|N} φ(d) k^{N/d} mod r (= N × #necklaces of length N in k colours, mod r)
  for small primes r chosen by the decoder and k a primitive root mod r;
* C: primitive periodic-point count of the cat map A = [[2,1],[1,1]] at period N, mod r: Σ_{d|N} μ(N/d)(L_{2d} − 2)
  mod r, r chosen by the decoder.

**Prediction P2:** decoders G, A, P, L split 100% (they are closed-form); K and C split ≥ 99% with ≤ 4n oracle
queries. Any failure is reported with its cause.

**E3. Point counts of X_0(N) over small fields leak quadratic symbols only.** n ∈ {32, 64}, 2000 instances each with
N ≡ 1 (mod 4). Oracle: #X_0(N)(F_l) = l + 1 − a_l(N) for primes l ≤ 13 (diagnostic, from (a) with the factors). By
the trace formula these values are functions of N and of the symbols χ_D(p), D = t² − 4l (proved in the results
note). Decoder: enumerate the 0/1 patterns of the unknown symbols, keep those consistent with every oracle value,
and output ε = (−1/p) if it is determined.

**Prediction P3:** (i) for every instance the true pattern is consistent (a correctness check); (ii) ε is
determined for ≥ 90% of instances (a measurement; we have no strong prior); (iii) no decoder output ever depends on
anything except symbols χ_D(p) of discriminants |D| ≤ 52, i.e. the oracle carries no information about p beyond
those ≈ 16 Legendre symbols. (iii) is a proved consequence of the formula; the experiment only checks the code.

## 3. Rules
* Exact integer / rational / F_P arithmetic only. Decoders (`operational`) never import factoring routines; oracles
  live in `src/geometry/harness.py` and are labelled diagnostic.
* Raw output: `experiments/raw/geometry_e{1,2,3}.jsonl`. Memory ≤ 2 GB; one job at a time.
* Deviations from this note are disclosed in the results note.
