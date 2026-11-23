# entropy: search the landscape of cheap N-only tests for persistent advantage on bits of p (2026-11-23)

**Motivation.** The user's question: can we explore tests by how much entropy about the factor they extract? The
right notion is computational entropy. Hardcore-bit amplification (ACGS 1988; Goldreich–Levin; Håstad–Näslund) means
one efficient test with advantage ≥ 1/poly(n), *persisting as n grows*, would factor. Predicted null: advantage decays
faster than any polynomial (Theorem PT, Theorem AC0-ε).

**Process.**
- Codex built and ran the pipeline (session 01a207c8-1b13-7cb1-8339-88e3d3bf936a).
- Gemini designed an independent grammar, targets and statistics protocol (agy 176da865):
  `experiments/entropy/entropy_spec.md`.
- Round 2 implemented Gemini's spec. Briefs: `codex/entropy_{common,codex,gemini,2}.md`. Outputs: `codex/entropy_*.out`.
- Code, pre-registrations, reports and tables: `experiments/entropy/round{1,2}/`. About 30 GB of raw data was not
  committed.

## Verdict: no lead. The only real signals decay exponentially, as N^{−1/4}.

**Round 1.**
- Search: 12,000 bit-valued tests on 1.6M semiprimes, n = 24–96.
- Discovery: all 482 survivors came from 24-bit numbers. That is a finite-population artifact: there are only 19,546
  distinct products, and every one also splits under a 256-step Fermat check. None replicated at 80 or 96 bits.
- Controls: the planted 1/n leak was detected; 1/n² was not.

**Gemini's specification (key contributions).**
- Statistical traps: universe exhaustion and prime reuse at small n; the power needed for 1/n² (about 2·10⁸ samples at
  n = 96); plug-in mutual-information bias; the log|adv| regression trap, which misclassifies exponential decay as
  polynomial 98% of the time; sorting bias from taking p = min(p, q); and size-normalization curvature.
- An exact size-free target, U_PNT(p|N) = 1 − ln(ln q/ln p)/ln((ln N − s_min)/s_min).
- A seven-family open-frontier grammar:
  - G1: deep and summed PRIM gates, modular-inverse gates, SLP-compressed powers.
  - G2: class-group NUCOMP powers beyond √N.
  - G3: multiplier continued fractions √(kN) and infrastructure jumps.
  - G4: X₀(N) / Frobenius proxies.
  - G5: small Coppersmith/GAR lattices.
  - G6: elliptic-curve ladders mod N.
  - G7: non-AC⁰ readouts.

**Round 2 (57 minutes).**
- Search: 20,000 truth-table-distinct tests across G1–G7 against 66 targets, corrected over 5.28M comparisons, on
  599.83M sampled products. The PNT encoder had a rounding bug, caught by audit, fixed and replayed byte-identically.
- **20 locked pairs were significant and split-clean at 64 bits. All 5 finalists passed Holm at 80 and 96 bits.** All
  five come from the first 64 continued-fraction steps of √(7N), i.e. three distinct comparisons:

  | Readout and target | Advantage at 64 / 80 / 96 bits (ppm) | Fitted β in 2^{−βn} | 2 log(L_poly/L_exp) |
  |---|---|---|---|
  | Even-index count parity → character of 2 | 450.8 / 29.3 / 1.89 | 0.246 | −6,416 |
  | Low-byte parity → character of 2 | 448.7 / 31.6 / 1.71 | 0.251 | −3,875 |
  | Jacobi readout → smaller factor mod 3 | 136.0 / 8.7 / 0.49 | 0.250 | −6,086 |

- **Interpretation.** These are genuine CF/genus leaks. Rank 1 fired on 959 of 224,971,804 clean rows at 96 bits
  (exact p = 2.6·10⁻⁴⁴). But the advantage decays as **2^{−n/4} = N^{−1/4}**, which is the square-form/infrastructure
  rate (SQUFOF/Lehman type), and the polynomial-decay hypothesis is decisively rejected. This is consistent with
  leakmap's conjecture H-GEN: continued-fraction statistics leak genus symbols, but only at rates that are exponentially
  small in n.
- **Power.** A top-up raised single-test power for a 1/n² leak at 96 bits to 82.3%. The plant measured +157 ppm, CI
  [92, 223], against 108.5 ppm planted. Power at the five-test Holm cutoff was 64.7%. The low-bit sentinel and shuffle
  controls were null. No PNT or Goldreich–Levin target locked.

**Bottom line.** The cheap-test landscape contains real but exponentially decaying leaks (the known archimedean/genus
family, at rate N^{−1/4}). No polynomially persistent advantage was found, so there is nothing to amplify.
