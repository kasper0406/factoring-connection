Explain, for a curious technically strong non-specialist (software engineer following this factoring project;
comfortable with modular arithmetic, CRT, groups, and "points of an elliptic curve mod p form a group"), the statement:

"The repo's Theorem EXT already recovers the factors from M = #E(Z/N) = (p+1−a_p)(q+1−a_q) for every good semiprime.
The only missing step is computing M without the factors. Schoof's algorithm, which counts points over a prime field,
has no known composite analogue."

Repo (read-only, do NOT edit): /Volumes/git/factoring-connection — notes/tracedata_results.md (Theorem EXT: read its
statement and proof; grep "EXT"), notes/ecount_results.md (Okamoto–Uchiyama; one count factors N),
notes/nativemod_results.md (why Edixhoven–Couveignes/Frobenius steps fail at composite N), notes/nonsmooth_results.md
(Prop FROB), notes/inp_results.md.
Write ~900–1300 words, plain language, short sections, minimal precise formulas:
1. E(Z/N) ≅ E(F_p) × E(F_q) by CRT; why its size is M; Hasse bound and the trace a_p.
2. The actual mechanism of Theorem EXT (faithful to the note: what is computed from M, the twist / auxiliary algebra if
   used, why a gcd reveals p, what "good" excludes), plus the simple intuition for why M carries p+q-type information.
3. How Schoof works over F_p (Frobenius π² − a_pπ + p = 0 on E[ℓ], division polynomials, x^p mod (ψ_ℓ, p), CRT over
   small ℓ) and exactly which step breaks over Z/N, and why "just run Schoof mod N" (using x^N) computes nothing useful.
4. Why this is the cleanest single missing step, and the honest caveat that computing M is exactly as hard as
   factoring (a reformulation, not a shortcut).
End with a 3-sentence summary. Final answer only (no file edits).
