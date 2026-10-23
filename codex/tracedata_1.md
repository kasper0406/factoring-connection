You are a research mathematician / computational number theorist. Blind round for agent `tracedata` (phase 25) of a
research project on classical integer factorization (repo: /Volumes/git/factoring-connection, read-only for you;
scratch files go in your current directory; python with sympy/gmpy2/python-flint/cypari2 is at
/Volumes/git/factoring-connection/.venv/bin/python). Rigour and honesty matter more than a positive result. Use web
search to check every literature claim; give exact citations (authors, title, venue, year, theorem/section number)
and mark anything you could not verify as UNVERIFIED. Label every claim KNOWN (with reference) / PROVED (proof given
here) / HEURISTIC / CONJECTURED / OPEN / FALSE. Negative literature findings must be stated as such.

Background (read first): notes/ecount_results.md in full (Theorems TW/TW*, the exceptional set X_E, §3 oracle
relations, §4 partial information). Optional: notes/nonsmooth_results.md §5.2 (Prop FROB), §5.5 (Prop RANK);
notes/factorbits_results.md §3.2 (Theorem PT); notes/leakmap_results.md §7 (Theorem PT-J).

Notation: N = pq, p < q primes (balanced unless stated), gcd(N, 6Δ_E) = 1. E/Q a FIXED elliptic curve (independent of
N), a_r = r + 1 − #E(F_r), a_N = a_p a_q (the Hecke eigenvalue of the newform at N), M = #E(Z/N) = (p+1−a_p)(q+1−a_q).

QUESTIONS (oracle questions: "X ⇒ FACT" means a (randomized / deterministic / ERH) polynomial-time algorithm that
factors N given N and the value X).

1. a_N alone. Given N and a_N = a_p·a_q for one fixed curve, can you factor N? (a) Count/describe the set of
   (p', q', a', b') with p'q' = N... — more usefully: for a FIXED N, how many factorisation-consistent explanations
   exist, and is the information-theoretic content of a_N about p nontrivial? (Hint to check, not to trust: the
   factorisation of N is unique, so the question is purely computational; |a_N| ≤ 4√N.) (b) Any reduction a_N ⇒ FACT
   (via factoring a_N itself? Galois-representation congruences a_p ≡ tr ρ_ℓ(Frob_p) mod ℓ? Coppersmith/lattices?
   several fixed curves E_1..E_k giving a_N(E_i)?) or an obstruction / evidence of hardness (e.g. a relativised or
   generic-model argument, or a reduction of a_N to/from something standard). (c) Bach–Charles (arXiv:0708.1192):
   state exactly what a(N²) supplies in their proof, and whether any single value a(N) or a(N^k) suffices.
2. The pair (a_p, a_q), ordered or unordered, or symmetric functions (a_p + a_q, a_p a_q). Does it factor N? What is the
   minimal symmetric trace data that provably factors N? Compare with M = N + 1 + s − (q+1)a_p − (p+1)a_q + a_pa_q
   (s = p + q): what does (a_p + a_q, a_p a_q, M) give, etc.? Is there a computable multiple of #E(F_p) or of the
   twist order from (N, a_p, a_q)?
3. Fixed-ℓ partial counts for a fixed curve: M mod ℓ (ℓ = 2, 3, 5, 7 fixed), and the number of roots mod N of a
   fixed polynomial (e.g. the cubic x³+ax+b, or ψ_3 of degree 4). Place each precisely: computable in poly time from
   N alone? Equivalent (both directions, poly-time) to fixed-base quadratic residuosity QR_N(Δ) for some Δ, or to a
   finite set of such bits, or to φ(N) mod ℓ? Factoring-hard? Note: the QR_N problem for a fixed base is not known
   to be factoring-hard. Give reductions in both directions where possible; where not, explain the obstruction
   (e.g. information in M mod ℓ that is not a function of quadratic symbols: non-abelian Frobenius data). In
   particular: is "number of roots of a fixed irreducible f mod N" for Galois group S_3 equivalent to QR_N(disc f)?
   For non-abelian f of higher degree (S_4 ψ_3, GL_2(F_ℓ)-type division fields), is the root count (or M mod ℓ)
   possibly HARDER than any finite set of quadratic-residuosity bits? Does the root count / M mod ℓ admit any
   poly-time algorithm (reciprocity-type) — compare with Theorem PT / PT-J and the abelian aggregate theorem.
4. X_E and determinism. ecount's TW* reduction (count ⇒ FACT) is rigorous off an explicit density-zero set X_E
   (M odd and both twist exponents divide M). (a) Give an everywhere randomized reduction for arbitrary fixed E using
   M (e.g. ℓ-adic descent at odd ℓ, points of E over Z/N[√d] or over cubic extension rings, the twist by a cubic
   character for j = 0, division-polynomial roots, or Weil-restriction/abelian-surface tricks), or prove that some
   natural class of methods fails on X_E. (b) A deterministic (possibly ERH) reduction from M to FACT for one fixed
   curve: small-height points, x-coordinates 1..(log N)^2, Bach-type bounds for least non-residues of f(x).
5. Literature: Bach–Charles; Kunihiro–Koyama EUROCRYPT'98; Martín–Morillo–Villar 2001; Okamoto–Uchiyama EUROCRYPT'98;
   Galbraith–McKee 2005; Dąbrowski–Pomykała–Shparlinski 2023; Edixhoven–Couveignes; also anything on "computing
   Fourier coefficients a(n) of modular forms at composite n is as hard as factoring", "trace of Frobenius modulo N",
   "point counting modulo composite", Schoof at composite modulus. Is any answer to Q1–Q4 already known?

Give a one-paragraph verdict at the top. Toy computations welcome (scratch dir only); report exact commands and numbers.
Write your full report as your final answer (saved to codex/tracedata_1.out).
