# wheelgap: short square-coset representatives on a primorial wheel vs a random-prime wheel (outside Brief 2, dir. 2) (2026-11-26)

**Process.** GPT-6.2-astra via codex (session 01a2176f-6e29-72b1-903a-efb366d6d884). Brief: `codex/wheelgap.md`.
Output: `codex/wheelgap_codex.out`. Prereg: `notes/wheelgap_prereg.md`. Code and CSVs: `experiments/wheelgap/`. The
brief's formulas were lost, so the problem was reconstructed.

**Verdict.** No sub-L[1/3] route. There is no primorial advantage. Codex gives the direction <1%.

- **Reconstruction (proved).**
  - The public sign vector is σ_W(N) = σ_W(p)σ_W(q). A rational a/b with max(a, b) ≤ H, H ≥ max(p, q), whose sign
    vector matches has Nab square once the wheel rejects every nonsquare in the box. Then a gcd splits N, with planted
    witness p/q.
  - Random-prime wheels (k = ⌈4 log₂ N⌉ primes of size exp((log N)³)) reject every nonsquare unconditionally,
    via the uniform PNT. Consecutive primes up to C(log N)² do so under GRH.
  - The "22-bit" clue is consistent with the first collision on the 168-bit 131#/2 wheel: 2,821,283 vs 3,617,267.
- **Literature (known).** This is prescribed-Legendre-pattern search: Lehmer and Shanks, pseudosquares,
  Lukes–Patterson–Williams, Sorenson. The scan ceilings are Graham–Ringrose for consecutive primes and Burgess
  W^{1/4+ε} in general. The brief's difference-table coefficient (log 2)/4 is the balanced half-wheel list size. It is
  the same constant as in Hittmeir's hyperbolic BSGS, and Bernstein's doubly focused enumeration gives (log 2)/2.
- **Counting obstruction (proved).** Fix an external base of k^{1+o(1)} primes and a height of 2^{O(k)}. Either wheel
  then reaches at most 2^{o(k)} of the 2^k sign vectors. This covers LLL, BKZ and exact CVP on that base equally for
  both wheel types. Larger or target-selected bases remain open. Thm J and Prop CV do not cover unrestricted wheel
  decoding.
- **Measured.** 42–168-bit wheels, 32 paired targets per size. For the exact shortest rational at 168 bits,
  log H/log W = 0.1093 (primorial) vs 0.1097 (random). The paired difference CI is [−0.0023, 0.0017]. LLL/Babai is
  far worse, at exponent ≈ 0.47. **Caveat:** the matched random wheel shared 29.9 of its 31 primes with the
  primorial, so it is a weak control. A 16-prime random control changes the constraint count, so it is not
  comparable.
- **Open.** A decoder that maps the primorial signs of a balanced pq to a factor-bearing rational of height O(√N) at
  cost L[α < 1/3].
