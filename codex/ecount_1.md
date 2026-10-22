You are a research mathematician / computational number theorist. Blind round for agent `ecount` (phase 23) of a
research project on classical integer factorization (repo: /Volumes/git/factoring-connection, read-only for you;
scratch files go in your current directory; python with sympy/gmpy2/python-flint is at
/Volumes/git/factoring-connection/.venv/bin/python). Rigour and honesty matter more than a positive result. Use web
search to check every literature claim; give exact citations (authors, title, venue, year, theorem/section number)
and mark anything you could not verify as UNVERIFIED. Label every claim KNOWN (with reference) / PROVED (proof given
here) / HEURISTIC / CONJECTURED / OPEN / FALSE. Negative literature findings must be stated as such.

Background (optional reading): notes/selfred_results.md §0 item "Most decisive open question" and §6 (a);
notes/nonsmooth_results.md §1 (equivalence table), §3, §5.2 (Prop FROB: producing the Frobenius operator on an étale
Z/N-algebra, e.g. Schoof's (Z/N)[x]/ψ_ℓ when (N/ℓ) = −1, factors N), §5.5 (Prop RANK); notes/factorbits_results.md
§3.2 (Theorem PT, "pure" generic-ring decision trees: information about the factors costs splitting).

Notation: N = pq, p < q primes of equal bit length, gcd(N, 6Δ_E) = 1. E/Q a FIXED elliptic curve y² = x³ + ax + b
(independent of N). #E(Z/N) := #E(F_p)·#E(F_q) (projective points of the Z/N-scheme; CRT). a_p = p + 1 − #E(F_p).

QUESTIONS

1. Dieulefait–Urroz, arXiv:1911.11004 ("Factorization and malleability of RSA moduli, and counting points on
   elliptic curves modulo N", or similar title; also published in Mathematics (MDPI) 2020?). Fetch it. State EXACTLY
   every theorem in the paper about #E(Z/N): hypotheses (single curve? curve plus twist? ERH? which definition of
   E(Z/N) — product of the local counts, or affine solutions plus one point at infinity?), conclusion
   (deterministic/randomized; does it need #E(Z/N) only, or the FACTORISATION of #E(Z/N)?), and the proof
   mechanism (Coppersmith? group-exponent/Miller-type argument?). Check each proof; report any gap.

2. One fixed curve. Is the map N ↦ #E(Z/N) (E fixed, ordinary, non-CM) factoring-equivalent?
   (a) Direction count ⇒ factoring. Given N and M = #E(Z/N) only (no factorisation of M), can you factor N in
       (randomized or deterministic, possibly ERH) polynomial time? Prove it or identify the obstruction precisely.
       Note: E(Z/N) has no efficiently samplable points in general without square roots mod N, unless you find a way
       around that. Characterise any exceptional set of (p, q) where your method fails and estimate its size; does
       it contain N for non-CM curves? What about CM curves?
   (b) Direction factoring ⇒ count: Schoof + CRT (KNOWN). Is anything cheaper than factoring known for computing
       #E(Z/N) for a fixed curve? Discuss Schoof's algorithm run over Z/N, Hecke operators T_N, modular symbols,
       Edixhoven–Couveignes, p-adic point counting, and the relation between #E(Z/N) and the modular coefficients:
       a_N = a_p·a_q (multiplicative), versus #E(Z/N) = (p+1−a_p)(q+1−a_q). Bach–Charles (arXiv:0708.1192):
       computing a(N), a(N²) of a fixed newform/τ factors N — exact statement? Is #E(Z/N) equivalent to the pair
       (a_p, a_q), to a_N, or strictly weaker/stronger (as oracles, given N)?
   (c) Lattice view. Given N and M = (p + B)(q + A) with unknown |A|, |B| ≤ 2N^{1/4}+1, is there a Coppersmith /
       Jochemsz–May / implicit-factoring style algorithm recovering p? What do the known bounds predict?

3. Partial information. For fixed E and small fixed prime ℓ: can #E(Z/N) mod ℓ, or the number of roots modulo N of
   the ℓ-division polynomial ψ_ℓ (or of x³ + ax + b), or the 2-torsion structure, or any other partial count be
   computed in polynomial time without factoring? What do these values reveal about p (e.g. via the Galois image
   mod ℓ), and does computing them reduce to or from known open bits (fixed-base quadratic residuosity, φ(N) mod ℓ)?
   Is it consistent with Theorem PT (pure trees learn nothing without splitting) and Prop RANK?

4. Give a one-paragraph verdict: what is known, what you proved, what remains OPEN.

Write your full report as your final answer (it is saved to codex/ecount_1.out). Toy computations are welcome
(scratch dir only); report exact commands and numbers.
