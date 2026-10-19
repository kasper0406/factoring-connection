"""detqs: deterministic quadratic sieve on almost all N (phase 8, polyimpl D1).

Modules
- sieve.py     OPERATIONAL: exact smoothness sieve of Q_t = (m+t)^2 - N, relations.
- linalg.py    OPERATIONAL: incremental exact F_2 elimination (Python-int bitsets).
- qs.py        OPERATIONAL: deterministic QS driver (sees only N and public config).
- localmodel.py MODEL: local-factor (small-prime) model of E[Z_N | N mod l], Dickman rho.
                Floating point is used here only for model predictions (fits), never
                for any measured quantity.
- harness.py   HARNESS: seeded instance generation (uses primality of p, q) and the
                pre-registered experiments E1-E3. Never imported by operational code.
"""
