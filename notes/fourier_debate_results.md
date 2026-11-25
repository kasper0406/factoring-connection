# fourier_debate: OpenAI "The logarithmic exponent of finite Fourier transforms" (2026-11-24) vs factoring

**Paper.** https://github.com/openai/math/blob/main/preprints/The-logarithmic-exponent-of-finite-Fourier-transforms-November-24-2026/paper.pdf
It shows that the linear circuit complexity of finite abelian Fourier transforms, and of formal polynomial
multiplication, is N(log N)^{1/2+o(1)} over every suitable field.
- Lower bound: a joint rank+support "rectangle" test.
- Upper bound: batched parity prediction (sign shifts give binomial concentration, a degree-√(m log m) polynomial
  gives the parity exactly, and exceptions at known locations are repaired) on prime tensor powers.
- Field reduction: Good/CRT axes, Rader, and shifted smooth primes p−1 via Bombieri–Vinogradov.

**Process.** A debate with GPT-6.2-astra via codex (session 01a21444-08c6-7c93-ad9f-ed48d84be183), two rounds. Briefs:
`codex/fourier_debate_{1,2}.md`. Outputs: `codex/fourier_debate_*_codex.out`. Experiment: `experiments/fourier_debate/`.

## Verdict: no sub-L[1/3] connection

- **Orchestrator's view.** Primes enter only as standard levers on lengths whose factorization the algorithm chooses
  itself. The direct payoffs (faster polynomial multiplication and multipoint evaluation for Pollard–Strassen,
  Bernstein batch smoothness trees, NFS linear algebra) are polylog, i.e. L^{o(1)}.
- **Codex round 1 ("qualified yes", 0.2%).** Steelman:
  - predict-and-repair as a branch selector for M2 quadratic CRT decoding;
  - RM(1,r) tied parity as a structure for the R hatch;
  - the rectangle test as a barrier tool;
  - Theorem T as the specification of what a Fourier repair would need;
  - choosing the factor base like the paper's prime axes.
- **Round 2: codex conceded.** Its answers to the four objections:
  - **R1.** The aggregate pair-count check is local CRT density plus box discrepancy, with no cheap way to get the
    discrepancy.
  - **R2.** The paper repairs at KNOWN locations, whereas M2 needs the locations as part of the answer.
  - **R3.** The binary circuit needs 2 to be invertible, but valuations are taken mod 2. The rectangle test is
    barrier-direction only.
  - **R4.** It has no concrete joint check from a factor base with smooth p−1.

## Experiment: the pair-count guide for M2 (codex; pre-registered, held out)

**Method.**
- Score boxes by Σ binom(ω_B(Q(u,v)), 2), i.e. pair prime-root lattice counts. Bisect a 64×64 box, ranking children
  by the score's excess over the local-density prediction.
- Data: 192 fresh forms per size at 64, 80 and 96 bits, near the one-witness threshold. Every relation was verified
  exactly; all 576 outcomes were audited with a full-box sieve.

**Results.**

| Bits | Guided nodes / relation | Unguided nodes / relation | Shared-score CPU / relation ÷ matched divisor sieve |
|---|---|---|---|
| 64 | 655 | 696 | 1.44× [1.30, 1.61] |
| 80 | 346 | 393 | 1.55× [1.40, 1.71] |
| 96 | 343 | 399 | 2.05× [1.95, 2.16] |

- The guide visits 6–14% fewer branches than unguided search, a real finite-box signal. Branch-saving growth is
  uncertain: the 96/64-bit ratio of the savings is 0.913 [0.823, 1.005].
- **Its CPU cost per relation is worse than the matched sieve, and gets worse with size**: the 96/64-bit cost ratio is
  1.43× [1.26, 1.60].
- **Post-hoc diagnostic.** On the 199 forms with exactly one witness, the score picked the witness's child 66% (upper
  four splits) and 90% (lower four) of the time. Replacing the witness's own binom(ω, 2) term with its local
  expectation dropped this to about 49–51%. **The guidance came from the scorer already counting the witness's own
  prime incidences**, i.e. from doing the sieve. The check is retired.
