# genus2rh pre-registration (agent `genus2rh`, phase 6, 2026-10-18)

Committed before any code for the experiment below was written or run.

**Scope.** The main deliverable is a proof (notes/genus2rh_results.md). No experiment can prove
or refute an asymptotic statement. The single experiment below is an exact correctness check of
the *algorithmic* half of the proof: the separation proposition (genus-2 analogue of Lenstra 1987,
Prop. 2.6) and the reduction lemma for Cantor's algorithm over Z/N. It says nothing about the
analytic inputs (LPP II Thm 1.1 needs p > 8100; Younis Thm 1.3 is asymptotic).

## Code (to be written after this file)
- `src/genus2rh/cantor.py`: polynomial arithmetic over Z/N that raises `NonUnit(d)` whenever a
  leading coefficient c with 1 < gcd(c, N) < N must be inspected; Cantor composition/reduction for
  y^2 = F(X), F monic quintic; scalar multiplication.
- `src/genus2rh/hecm.py` (**operational**: sees only N, public parameters v, w and a seed): sample
  (a, b, f1..f5) uniformly in (Z/N)^7, set f0 = b^2 - sum f_i a^i, check gcd(f5, N) and
  gcd(disc, N), pass to the monic model, D = (X - A, B); for primes r <= w in increasing order,
  each repeated e(r) times (e(r) = max m with r^m <= (sqrt v + 1)^4), replace D by [r]D.
  Return the first proper divisor met, else None.
- `src/genus2rh/diagnostic.py` (**diagnostic**, factor-aware): #J(F_p) by exact point counts over
  F_p and F_{p^2}; ord(D_p) by factoring #J(F_p) (sympy); the conditions below.
- `src/genus2rh/experiment_g1.py`; raw output `experiments/raw/genus2rh_g1.jsonl`.

## Experiment G1 (exact implication check)
- Instances: 40 pairs p < q of primes drawn uniformly (seed 20261018) from [101, 251]; N = pq.
  For each N, 30 curve tuples (seeded per N). Parameters v = 251, w = 50.
- For each run and each ordered pair (r, r') in {(p, q), (q, p)}, compute the conditions
  (ii) f mod r is a squarefree quintic; (iii) #J(F_r) is w-smooth; (iv) f mod r' is a
  squarefree quintic; (v) the largest prime l dividing ord(D_r) does not divide #J(F_r').
  Call the run COVERED if (ii)-(v) hold for at least one ordered pair.
- **P1 (prediction, a proved implication):** every COVERED run returns a proper divisor of N.
  Expected violations: 0. Any violation refutes Proposition A as implemented, and I report it.
- **P2 (prediction, the reduction lemma):** for every run that returns None, the final Mumford
  pair of [k]D over Z/N, reduced mod p and mod q, equals the result of the same computation done
  over F_p and over F_q. Expected mismatches: 0.
- **Descriptive (no prediction):** the COVERED rate; the success rate among non-COVERED runs;
  setup-stage divisors (gcd with f5 or with the discriminant: ACCIDENTAL, null rate about
  2/p + 2/q... per run) vs main-loop divisors; and for main-loop divisors, the fraction for which
  neither ord(D_p) nor ord(D_q) divides k (accidental intermediate coincidences).
- Also checked in the unit tests (not part of G1): #J from point counts annihilates random
  divisors under F_p Cantor arithmetic, and the group law is associative and commutative on
  random triples.
