# prim pre-registration (phase 25)

Agent `prim`, 2026-10-23. Committed before any experiment code exists and before any experiment is run. Code will live
in `src/prim/`, raw data in `experiments/raw/prim_*.jsonl(.gz)`, summaries in `experiments/analysis/prim_*.json`.

## Theory under test (to be proved or labelled in `notes/prim_results.md`)

Notation: N = pq, p, q primes in I = (X, 2X] (b-bit primes with the top two bits set), n_p = p − 1. A node v is a
computed element; its answer at p is [v ≡ 0 mod p]. Rank-1 algebra A = Z unless stated.

- **Δ-form / uniform-r domination (U).** For an affine menu exponent E(N) = αN + β, at p the gate y^{E(N)} equals
  y^{αr+β} in F_p with r = N mod (p − 1) = q mod (p − 1). Every residue r mod (p − 1) is hit by at most 2 primes q of
  I, so for every set Z ⊆ Z/(p − 1): #{q : r_q ∈ Z} ≤ 2|Z|. Hence the deviation frequency of a node at p is controlled
  by the *uniform-r zero set* Z_v(p) = {r ∈ Z/(p−1) : node v vanishes at p when every gate exponent is αr + β}, a pure
  function of p.
- **KS (Kelley structure, one-layer summed gates).** If at p the node is f(g^r) with f a t-sparse Laurent polynomial
  (t ≥ 2), then Z(f) = V ∪ R where V is a union of cosets of subgroups of order ≥ n^{1/e} on each of which f vanishes
  identically (order alignment, (M)) and |R| ≤ 2n^{1 − 1/(e(1 + ⌊log₂ t⌋)(t − 1))} (to be proved; the exact constant in
  the results file is authoritative if it differs; P4 checks whichever version is proved).
- **H-RAND (heuristic under test).** For nested gates with deterministic pure bases, the zero set of
  A(r)·u(r)^{e(r)} − B(r) (u(r) the shifted inner value, e(r) = αr + β) behaves like the *fresh-random model* in which
  u(r) is replaced by an independent uniform element of F_p^*: Pr[zero at r] = m(r) := [B/A ∈ (F_p^*)^{d}]·d/n,
  d = gcd(e(r), n). The fresh-random model is (M) by a Weil argument (to be proved for depth 2). H-RAND fails if
  nesting manufactures structure beyond (M)/(A).

## E1: instrumented PRIM programs on small prime pairs

Populations (unordered pairs p < q, both b-bit with the top two bits set; iid law = uniform distinct pairs):
- b = 12 and b = 14: all pairs. b = 16: 200 000 pairs sampled uniformly (seed `primE1|16`); full uniform-r zero sets
  for every prime of the b = 16 range.

Nodes (each a zero test of the displayed element; affine menu unless marked floor):
- Nested (deterministic bases): D2a_w: (2^N + 1)^N − w, w ∈ {1, 2, 3}; D2b: (2^{N−1} + 1)^{N−1} − 1;
  D2c: (3^{N+1} + 2)^{N−1} − 4; D2d: 2^N·(2^N + 3)^N − 5 (retained product); D2e: (2^N + 3^N)^{N+1} − 1 (summed then
  gated); D3: ((2^N + 1)^N + 1)^N − 1 (three layers); D2f (floor): (3^{(N−1)/2} + 2)^{(N+1)/2} − 1 (q-world only);
  L1 (A = Z[x]/(x² − 3x + 1), α = x, element test): (α^{N+1} + 2)^{N−1} − 1 (q-world only).
- Summed (one layer): S1: 2^N + 3^N − 5^N; S2: 2^{N−1} + 3^{N−1} − 2; S3: 2^N + 3^{N+1} + 5^{N−1} + 7;
  S4: 2^N − 2·3^N + 3·5^N − 4·7^N + 5·11^N − 6·13^N + 7; S5: (2^N + 1)(3^N + 1) − 2.
- Controls: C_F: 2^{N−1} − 1 (Fermat, M); C_A_k: 2^N − 2^k, k = 1..16 (archimedean pins, A);
  C_P: (2^N + 1)^M − 1 with M = lcm(1..64) (gate followed by a pure smooth power, M); C_PM1: 2^M − 1 (pure, PT).
- Batteries (non-adaptive; transcript = vector of answers): NEST = {D2*, D3, D2f}; SUM = {S1..S5}; CTRL = controls;
  ALL = union; LUC = {L1}.

Reference κ(p) for PT-X: per-node majority answer over the partners of p in the population (a function of p alone).

Deviation classes for a zero at p of a node whose majority is "nonzero" (r = q mod (p − 1)):
- COSET (M): r lies in a coset r₀ + mZ/(p−1) (m | p − 1) of size (p−1)/m ≥ L := ⌈(p−1)^{1/3}⌉ contained in Z_v(p);
- PIN (A): otherwise, if q − p ≤ 2^{b/2} or r ≤ 64;
- OTHER: otherwise.
SMALL-flag at p: some gate base of the node has ord_p ≤ 2^{b/2}.

Targets: ε = (−1|p) on N ≡ 1 (mod 4); [p ≡ 1 (mod 3)] on N ≡ 1 (mod 3); [q − p < X/8].

Registered predictions:
- **P1 (diagonal lemma):** for every pair and node, the run mod N splits iff the answers at p and q differ, and
  otherwise returns the common answer; 0 mismatches (all b).
- **P2 (PT-X):** with κ as above and Bayes-optimal leaf labels computed exactly on the population,
  Pr[correct or split] ≤ B + 2s + 2Pr[D ∧ no split] in every (b, battery, target) cell.
- **P3 (U):** for every p and node, #{partners q : node deviates at p} ≤ 2·|Z_dev,v(p)| (affine nodes, A = Z);
  0 violations.
- **P4 (KS):** for S1–S5 and C_F at every p: every zero lies in a vanishing coset of size ≥ (p−1)^{1/e} or in the
  residual, and the residual obeys the proved KS bound; 0 violations.
- **P5 (H-RAND, full-r):** for each nested affine node (D2a_w, D2b–D2e, D3), over primes that are not SMALL-flagged:
  Σ_p |Z_v(p)| / Σ_p Σ_r m_p(r) ∈ [0.75, 1.33] at b = 14 and b = 16; and no non-flagged p has a vanishing coset of size
  ≥ L in Z_v(p).
- **P6 (H-RAND, q-world):** for each nested node (affine, A = Z), the number of pairs with a zero at the smaller or
  larger prime, divided by the summed model probability, lies in [0.6, 1.6] (b = 14, 16). Two-sided zeros (N divides the
  node, non-split deviations) summed over the nested nodes: observed ≤ 3·(Σ model products) + 5.
- **P7 (summed nodes, residual):** for S1, S3, S4, S5 the residual |R_p| (zeros outside vanishing cosets of size ≥ L)
  has mean over p in [0.3, 3] and max ≤ 20 (Poisson(≈1) heuristic for zeros of a random-looking function).
- **P8 (Jacobi):** on N ≡ 1 (mod 4) pairs, the dyadic-robust z of the ε-correlation of J = (c|N) is < 4 in absolute
  value for every J below, at b = 14 and b = 16; the positive control J_CRT (c built by CRT to be a QR mod every
  ε = +1 prime of I and a QR mod a fixed half of the ε = −1 primes) has z > 10 at b = 12 and b = 14.
  J nodes: J1 = (2^N + 1 | N); J2 = ((2^N + 1)^N + 3 | N); J3 = (2^{(N−1)/2} + 2 | N); J4 = (3^N + 2^{N+1} + 1 | N);
  J5 = (2^{N−1} + 3 | N); J6 = (V_N + 1 | N) with V_N = α^N + α^{−N} in Z[x]/(x² − 3x + 1);
  J7_k = (F_k mod N | N), F_k = 2^{2^k} + 1, k = b, …, b + 3 (pure compressed); J8 = (2^M − 3 | N), M = lcm(1..64).
- Every split found is classified STRUCTURAL (answers differ at p and q) or ACCIDENTAL; ACCIDENTAL expected 0.

Analysis decisions fixed now: dyadic-robust variance = observations dependent iff they share a prime (either role),
as in leakmap v2. Pass/fail is reported for every registered prediction; nothing registered will be dropped.
