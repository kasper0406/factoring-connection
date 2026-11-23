# entropy round 2: implement Gemini's spec, and fix power

Gemini's independent spec: /Volumes/git/factoring-connection/codex/entropy_spec_gemini.md (and codex/entropy_gemini.out;
verification script codex/scratch/ent_g/verify_entropy_spec.py). Key points to adopt:
- Statistical traps: small-n universe exhaustion / prime reuse (avoid n ≤ 32 for discovery); n=96 power (1/n² needs
  ~2·10^8 samples per locked test); plug-in MI bias; the log|adv| OLS trap (use signed advantage and a likelihood-ratio
  test of polynomial vs exponential decay on signed estimates); sorting bias of p = min(p,q); size-normalization
  curvature — use Gemini's exact PNT probability-integral-transform target U_PNT(p|N) for size-free top/ratio bits.
- Grammar G1–G7 (open frontier): deep/summed PRIM gates, modular-inverse gates, SLP-compressed powers, extension
  coordinates; class-group NUCOMP at exponents > √N for disc −4kN; multiplier CFs √(kN) and infrastructure jumps;
  X_0(N)/Frobenius proxies; small Coppersmith/GAR lattices; elliptic-curve ladders mod N ([N]P = [(a_p−1)q]P coupling);
  non-AC⁰ readouts (dense parities, base-b digit sums, automata).
Plan (pre-register before running; this session has your round-1 pipeline):
1. Implement as much of G1–G7 as runs ≤ ~50 µs/eval at n=64 (C++ where needed); dedupe by truth-table hashing.
2. Discovery at n ∈ {40, 48, 56} with fresh, non-repeating semiprimes and Gemini's targets; strong multiple-testing
   control; lock the top ≤ 20 survivors (plus planted 1/n and 1/n² controls).
3. Confirmation with large samples on fresh data: m ≈ 4·10^7 at n = 64 for all locked tests; m ≈ 2·10^8 at n = 80 or 96
   for the top ≤ 5 if compute allows. Report signed advantage with exact CIs per (test, n), the polynomial-vs-
   exponential LRT, planted-control power actually achieved, and accidental-split flags.
Budget ≤ 3 hours. Files in cwd; final summary ≤ 1500 words.
