You are a careful research mathematician and literature expert in computational number theory. Use web search. Repo context (read-only): /Volumes/git/factoring-connection, especially notes/dcltheorem_results.md (Theorem D′).

QUESTION (prior art only). Has anyone published or posted (journal, conference, arXiv, thesis, ePrint, blog, MathOverflow) a result of this form:

  "a DETERMINISTIC algorithm that always outputs the complete factorisation of N, whose cost AVERAGED over balanced semiprimes N = pq (p < q ≤ 2p, or p ~ q, or the RSA-type distribution of uniformly random n/2-bit primes) is N^{c} with c < 1/5"

i.e. an average-case (over inputs, not over coins) deterministic bound below the Harvey / Harvey–Hittmeir worst-case exponent 1/5 (Math. Comp. 90 (2021), 91 (2022))? Also accept: "for all but a small fraction of balanced semiprimes, deterministic time N^{c}, c < 1/5" (Theorem D′'s form), and average-case statements for Lehman / Fermat / McKee / McKee–Pinch / Hittmeir / SQUFOF / Pollard rho with fixed polynomial / Pollard p−1 / class-group (Shanks, Schnorr–Lenstra) / deterministic ECM with fixed curves, if they are PROVED (not heuristic).

Specifically check:
1. McKee–Pinch 1996 "Old and new deterministic factoring algorithms" (ANTS II): any average-case PROVED statement?
2. McKee 1999 "Speeding Fermat's factoring method" (Math. Comp.): heuristic N^{1/4} — any proved average version?
3. Hittmeir (2018, 2021, 2023), Harvey (2021), Harvey–Hittmeir (2022, 2026 order-selection arXiv 2601.11131), Gao–Feng–Hu–Pan arXiv 2512.19076, Stănică–Mulder–Hittmeir arXiv 2606.13018: any remark on average-case or early-exit cost?
4. Carella arXiv 1308.2891 (claims deterministic 1/6 worst case).
5. Any "average-case analysis of deterministic factoring" (e.g. Knuth–Trabb Pardo style analyses of trial division, Lehman average, "Fermat's method is fast on average"?).
6. Pollard p−1 with fixed base and deterministic stage 2: proved average over primes p (smoothness of p−1) — what exponent does it give on balanced semiprimes? (We believe it cannot beat 1/5 on average because the failure probability is a positive constant, not a power saving; confirm or refute.)
7. Anything on "deterministic factoring for almost all integers / almost all semiprimes" in sub-N^{1/5} time.

Output: a curated list with exact citations and URLs; for each, one line on whether it states a PROVED average-case or almost-all deterministic exponent < 1/5. Label each item KNOWN / NOT FOUND / UNCERTAIN. End with a one-paragraph verdict: is "deterministic complete factoring with average cost N^{1/6+o(1)} (even N^{1/6} polylog) over balanced semiprimes" stated anywhere you could find? Do not claim novelty; say "not found by us" at most. Keep it under 1500 words.
