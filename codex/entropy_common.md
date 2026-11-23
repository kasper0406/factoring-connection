# entropy: search the landscape of cheap tests on N for ANY persistent advantage on bits of p

Motivation (user): explore tests by how much computational entropy about the factor they extract. Theory: given N, p
has zero Shannon entropy but (if factoring is hard) full computational entropy. Hardcore-bit results (Alexi–Chor–
Goldreich–Schnorr 1988; Goldreich–Levin; Håstad–Näslund) mean a single efficient test predicting a bit / parity of p
with advantage ≥ 1/poly(n) PERSISTING as n grows could be amplified (list decoding) to factoring. So the experiment
looks not for big leaks but for advantage that decays only polynomially in n = log₂N. Predicted null: advantage
decays faster than any polynomial (or is exactly 0 for tests that are functions of N used only as a modulus — Theorem
PT, notes/factorbits_results.md; Theorem AC0-ε, notes/leakmap_results.md §6.1).

Known trivial leaks to control: (i) SIZE — N's magnitude constrains top bits of p and q (leakmap: AC⁰ predicts
min(p,q) > 3X/2 with constant advantage); condition/residualize on N's top bits or define targets that are size-free;
(ii) LOW BITS — p mod 2^k is uniform given N mod 2^k (pq ≡ N), so low-bit tests should show 0 advantage — a good
sanity check; (iii) ACCIDENTAL SPLITS at toy sizes — tests like x^E mod N, Lucas sequences, CF of √N can factor N
outright when p−1 / p+1 is smooth or p ≈ q (Lehman/Fermat); every instance where any test's intermediate values share
a factor with N (gcd > 1) or where Fermat/Lehman/p±1 at the test's budget would split must be flagged and the
advantage reported both including and excluding them (leakmap found archimedean leaks only via Lehman-type splits).
Targets (p = smaller factor, balanced semiprimes, both primes uniform in (X, 2X]): individual bits of p at several
positions (top after removing size info, middle, low), random GL parities ⟨r, p⟩, ε(p) = χ₄(p), Legendre symbols
(ℓ|p) for small fixed ℓ, p mod small m, and the "ratio" bits of p/q.
Repo: /Volumes/git/factoring-connection (do NOT edit; scratch only in your cwd). Python: repo .venv has gmpy2, cypari2,
numpy, sympy, flint; make your own uv venv in cwd for anything else (UV_CACHE_DIR in cwd).
