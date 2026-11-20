# guide: is there ANY predictive signal for smoothness beyond size and local congruence densities?

Motivation (user's idea): "plant seeds and tame their growth, tracking densities with another system as guidance".
The α ≥ 1/3 model theorem (notes/lowerbound_results.md) already allows arbitrary adaptive guidance UNDER assumption H
(conditional Dickman: given everything computed so far, a fresh probe of size s is B-smooth with probability
≤ ρ(s/log B)^{1−o(1)}). Guidance can only help if it exploits information H assumes away. Local congruence
information (which classes are divisible by small primes) is the divisor trick — proved to trade one-for-one
(Theorem J, weighted J; notes/hatch_results.md, notes/chase_results.md). Size is the lattice-geometry question
(KH floor). The remaining possibility is a GLOBAL correlation: smoothness predictable from features beyond size and
local divisibility. This experiment tests that directly. A positive result would be the first evidence against H.

Core design principle (no leakage): fix a small-prime cutoff y. For each candidate (an NFS (a,b) pair, or a
sublattice point in lattice-sieve coordinates), the TARGET is whether the y-rough cofactor (after removing all primes
≤ y on both sides) is B-smooth on both sides (and, separately, one side). The BASELINE predictor uses only (i) the
log-sizes of the y-rough cofactors and (ii) exact local information available from a sieve up to y. TEST features
must not reveal divisibility of THIS candidate by any prime in (y, B] (note: neighbours' divisibility by p in (y,B]
leaks via periodicity a ≡ a+p — that is local sieve information and must be excluded or residualized). Candidate test
features: positions (a, b, gcd structure), lattice coset labels mod products of primes > B, outcomes of OTHER
candidates far away / in other cosets, search history, smoothness of related algebraic elements (conjugates, a±1,
products), the big-prime part's digit/bit statistics of other candidates, and anything a learned model can extract
from raw (a,b) bits. Evaluate on held-out N and held-out regions; pre-register metrics (log-loss / AUC improvement over
baseline with CIs) and the decision threshold.

Repo / tools: /Volumes/git/factoring-connection (do NOT edit; scratch only in your cwd). CADO-NFS is built at
/Users/knielsen/src/cado-nfs/build/kmini (sieve/las, polyselect) — notes/cadobench_results.md describes earlier use.
Python: the repo .venv has gmpy2, cypari2, sympy, mpmath, flint, numpy (no sklearn/lightgbm); create your own venv in
your cwd with uv if you need ML libs (set UV_CACHE_DIR inside your cwd).
