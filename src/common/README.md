# src/common — shared infrastructure (owner: infraA)

Import as `from src.common.ring import ...` with the project root on `sys.path`
(run scripts as `.venv/bin/python -m ...` from the root, or `pytest` from the root).
The API below is intended to stay stable; additions only.

## metrics.py
- `Metrics(**fields)`: `inc(key, n=1)` additive counters; `hi(key, v)` running maxima;
  `set(key, v)` free fields; `start(name)/stop(name)` timers; `to_record(**extra)` -> dict
  with `counts`, `maxima`, `wall_s`. `NullMetrics()` ignores everything.
- `append_jsonl(path, record)`, `read_jsonl(path)`.
- Standard keys: `rn_mul`, `rn_add`, `gcd`, `inv`, `q_mul`, `q_sqr`, `q_mul_dimsum`,
  maxima `quot_dim`, `ext_degree`, `max_operand_bits`, `max_degree`, `stored_coeffs`,
  `circuit_nodes`.

## ring.py (operational-safe: never uses factors)
- `NonUnit(divisor, modulus, value)`; `.is_proper` iff `1 < divisor < N`.
- `inv_mod(a, N, m=None)`: gcd-checked inverse, raises `NonUnit`.
- `Zmod(N, m)`: R_N; elements are ints in `[0,N)`. Methods `add sub neg mul inv pow
  is_zero eq from_int`, attributes `zero one N`.
- `PolyQuotient(base, f, m)`: `base[X]/(f)`, f monic, list low->high of base elements.
  Elements are tuples of length `deg f`. Methods `add sub neg mul sqr pow scale reduce
  from_coeffs from_base from_int gen coeffs is_zero eq dim total_dim`.
  Towers work (`base` may itself be a `PolyQuotient`). Pure-Python reference.
  - `truncated(base, K)` = `base[X]/(X^K)`
  - `binomial_quotient(base, r, c)` = `base[X]/(X^r - c)`
  - `ext_ring(base, g)` = `base[t]/(g)`
- `FlintQuotient(N, f, m)`: R_N[X]/(f) (f ints, monic) via python-flint `nmod_poly`
  (N < 2^62) or `fmpz_mod_poly`. Same element interface (elements are flint polys,
  use `.coeffs(a)` to get a list of ints of length d). Only base = R_N.
  Counts `q_mul/q_sqr/q_mul_dimsum` only (flint's inner R_N ops are not counted).
- `power_of_linear(ring, e, a=None, b=None)`: `(a+bX)^e` by repeated squaring
  (default `(1+X)^e`).
- `flint_binomial_power(N, r, c, e)`: coeffs of `(1+X)^e mod (N, X^r - c)`.
- `flint_audit(bits=4096)`: empirical check that flint does not factor the modulus.

**Flint audit (verified 2026-10-17, python-flint 0.9.0).** We only call
multiplication, `mul_low`, `%` by a monic polynomial, shifts and truncation.
None of these needs an inverse other than of the leading coefficient 1.
(a) Results agree with `PolyQuotient` on random inputs and on `(1+X)^e` for
composite N of 10, 40 and 71 bits (tests). (b) A context over a 4096-bit RSA-shaped modulus is
created in ~10 ms, and `(1+X)^N mod (X^64-3)` is computed in ~3 s, so no
factorisation can be happening. (c) `fmpz_mod_ctx.is_prime()` is evaluated only
when called; operational code must never call it, nor `factor`, `gcd`, `xgcd`,
`inverse_mod`, `roots` on composite-modulus polynomials. Those routines assume
a prime modulus and may raise or behave like a hidden factor search.

## instances.py (HARNESS: knows the factors)
`Instance(N, family, params, seed, labels)`: candidates see `inst.N` / `inst.public()`;
`inst.labels` (`factors`, `p`, `distinct`, ...) is harness-only and hidden from `repr`.
Generators (deterministic in `(family, params, seed)` via SHA-512 string seeding):
`balanced_semiprime(bits, seed)` (bits in `BALANCED_BITS = 16,24,32,48,64,96,128`),
`unbalanced_semiprime(bits, small_bits, seed)`, `prime_power(bits, k, seed)`,
`multi_prime(bits, r, seed)` (r = 3, 4), `carmichael_small(index)`
(`CARMICHAEL_SMALL` list), `carmichael_chernick(kbits, seed)`, `prime_instance(bits, seed)`,
diagnostic `diag_near_equal(bits, seed)`, `diag_smooth_pm1(bits, seed, B=50)`.
`make(family, ...)`, `batch(family, count, ..., seed0=0)`. Also `is_prime`, `random_prime`.
The full sampling procedure is in the module docstring.

## baselines.py (operational)
`trial_division(N, m=None, bound=None)`, `pollard_rho_brent(N, seed=0, m=None,
max_iters=2^26)`: return a proper divisor or None; count `trial_div`, `rn_mul`, `gcd`.

## Elsewhere
`src/seed_checks.py`: BRIEF §9 smoke test verbatim + `check_extended(limit)`,
`kummer_support(N)`, `lucas_support_squarefree_part(N, p)`, `dominated`, `carries`
(reference code, uses trial division).
