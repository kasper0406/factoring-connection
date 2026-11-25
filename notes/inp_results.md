# inp: "suppose factoring ∈ P": where must a classical poly-time algorithm live? (2026-11-25)

**Context.** The user reported rumours (unverified; we have no information about any such result) that factoring is
in P. This is a working-hypothesis analysis: a debate with GPT-6.2-astra via codex (session
01a21498-3e5f-7d92-a799-a290d1f75bbb), two rounds. Briefs: `codex/inp_{1,2}.md`. Outputs: `codex/inp_*_codex.out`.
Pilot: `experiments/inp/`.

## Key structural point (codex)

A poly-time factorer need not find a "fourth leak family". It can produce a residue divisible by exactly one of p, q
and finish with a gcd. That is consistent with Theorem PT. Only *relation* methods are forced to break one of the
model theorem's assumptions: H, the size profile, or R.

## Ranking after Claude's pushback (conditional on factoring ∈ P; subjective)

| Weight | Location | Note |
|---|---|---|
| 24% | **Native fixed-curve point count over Z/N** | Recovery from M = (p+1−a_p)(q+1−a_q) already exists for every good semiprime (tracedata Theorem EXT). The single missing step is computing M mod ℓ at enough small ℓ at composite N; Schoof's prime-field Frobenius step has no composite analogue. **Codex's top route.** |
| 20% | Unidentified mechanism | The catalogue has no theorem covering arbitrary compositions. |
| 16% | Adaptive factor-interval refinement | A computable rule must keep the true factor among poly-many intervals down to width N^{1/4}, where Coppersmith finishes. Downgraded from 32% "mixed modular/integer", which is vacuous since integer representatives are computationally complete. |
| 13% | Exact class numbers / regulators | Exact h(−4N), or a certified multiple of hR, yields ambiguous forms (genusleak). The barrier is the nativeana precision wall. |
| 12% | Constructive relations | Must break H, the profile or R. |
| 8% | Fast modular factorial / dense arithmetic circuit | A gcd finish exists. The τ-conjecture argues against short integer circuits. |
| 4% | Structured groups | Random-input words in SL₂(F_p) are killed by Larsen–Shalev (≤ \|G\|^{−c}) and, for N-dependent poly-length words, by an O(ℓ/p) zero count. Survivors: N-dependent structured inputs or exponential-expanded-length programs, with no mechanism known. |
| 2% | Classical period finding | Theorem T. |
| 1% | A single non-CM coefficient a_f(N) | Two missing steps. The [N+1−A]P probe was retracted (no mechanism: a_N = a_p a_q is not a group order). Bach–Charles needs a_f(N²) as well. |

## Pilot run (codex): dual-Coppersmith sibling-interval score

**Method.** Every polynomial in (N, a+x)^m evaluates at p−a to a multiple of p^m. So the dual of the scaled lattice
contains p^{−m}(1, t, t², …) with t = (p−a)/X, and |t| ≤ 1 for the correct child interval. Use this to score which
child interval holds p.

**Result.**
- At radius about N^{0.12}: 216/216 correct, but that regime is already Coppersmith-easy.
- At N^{0.28}: 111/216. At N^{0.34}: 102/216. Chance is 108.
- **Zero factors in 432 above-threshold trials.**

The score is retired. The interval route would need a different invariant.

## Status

No sub-L[1/3] lead. The most concrete open target remains a native algorithm for #E(Z/N) of one fixed curve.
