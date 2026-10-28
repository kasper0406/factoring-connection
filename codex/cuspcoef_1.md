You are a research mathematician / computational number theorist. Blind round for agent `cuspcoef` (phase 30) of a
research project on classical integer factorization (repo: /Volumes/git/factoring-connection, read-only for you;
your cwd is a writable scratch dir; python with sympy/gmpy2/python-flint/cypari2 (PARI) is at
/Volumes/git/factoring-connection/.venv/bin/python). Rigour and honesty matter more than a positive result. Use web
search to check every literature claim; give exact citations (authors, title, venue, year, theorem/section number)
and mark anything you could not verify as UNVERIFIED. Label every claim KNOWN (with reference) / PROVED (proof given
here) / MODEL (proved inside a stated probabilistic model) / HEURISTIC / CONJECTURED / OPEN / FALSE. Keep memory use
under 2 GB and leave no background processes.

Optional background (read if useful): notes/nativemod_results.md (§0, §1.3, §3.2–3.3, §4, §5),
notes/tracedata_results.md (§0, §2.1, §2.5, §2.6 Theorem SC, §2.8), notes/ecount_results.md (§3).

SETTING. N = pq, p < q unknown odd primes, balanced (p, q in [X, 2X]), coprime to the level. f is a FIXED non-CM
newform with rational coefficients; main examples: the weight-2 forms of 37a1 (y^2 + y = x^3 - x) and 11a1, and
Delta (weight 12). a_f(N) = a_f(p) a_f(q). Equivalent forms of the datum: a_E(N) = a_p a_q = sum_{x mod N} ((x^3+ax+b)/N)
(Jacobi symbol); at level 11, a_{11a1}(N) = (r_{L1}(N) - r_{L2}(N))/4 for two quaternary lattices of one genus.

KNOWN: Bach–Charles (Contemp. Math. 463, 2008): a(N) and a(N^2) together factor a density-one set of RSA moduli.
Edixhoven–Couveignes: a(p) for prime p in poly(log p). Whether a(N) ALONE factors N is open; no native algorithm
(one that never factors N) faster than ~N is known to us. The problem is in ZPP^FACT.

QUESTIONS (answer each; negative answers with precise reasons are valuable).
Q1 (reduction). Is there a (randomized or deterministic, worst-case or almost-all) polynomial-time reduction from
  factoring N to (a) one value a_f(N) for one fixed f; (b) the values a_{f_i}(N) for O(1) fixed non-CM f_i;
  (c) for K = poly(log N) fixed non-CM forms (e.g. the first K curves y^2 = x^3 + Ax + B in a fixed enumeration)?
  Count information against constraints. Consider: Coppersmith/lattice attacks using |a_p| <= 2 sqrt p and
  a_N = a_p a_q (and several curves; and weight 2 together with weight 12, i.e. a_E(N) with tau(N)); Galois
  congruences (a_p mod l = tr Frob_p) plus Chebotarev; the integer factorisation of a_N; anything else. If (c) works
  only after factoring the integers a_N(E_i) (~ n/2 bits each), is there a way around that?
Q2 (native algorithm). Best native cost for a_f(N) at composite unfactored N: modular symbols (Merel/Heilbronn sets
  for T_N), theta differences (4-dim lattice-point counting), the complete Jacobi-symbol sum above, shifted
  convolution sums (r_{L1}(N) = sum_m r_B(m) r_B(N - m) with B = x^2 + xy + 3y^2 of class number one), analytic routes.
  Is anything below ~N^{1-delta}, let alone below L_N[1/3], known or plausible? Give the cost of each route.
Q3 (barrier). If neither, state and PROVE the sharpest model barrier you can: e.g. "in model M, a_N carries at most
  c bits about p", or "any reduction of type T from a_N to factoring yields a factoring algorithm". Be explicit about
  what the model allows the algorithm to compute from N (group law on E(Z/N), torsion algebras (Z/N)[x]/(psi_l),
  Jacobi symbols, ...) and which known reductions (Bach–Charles; two CM curves y^2 = x^3 - x, y^2 = x^3 - Dx factor
  N from a_N alone when p, q = 1 mod 4) the model correctly does or does not exclude.
Q4. Literature: any prior statement (positive or negative) on factoring from a single Hecke eigenvalue / Fourier
  coefficient at a composite index, beyond Bach–Charles, Charles (JCSS 2007, lattice counting #P), Kunihiro–Koyama,
  Dieulefait–Jiménez Urroz, Edixhoven–Couveignes Rem 2.2.4, A. Chow's thesis?

Write a structured answer (<= 2500 words) with the strongest results first. Small exact PARI experiments are welcome
(keep them under a few minutes each); report exact numbers and commands.
