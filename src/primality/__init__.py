"""Primality-test leakage study (owner: agent `primality`, Phase 2, 2026-10-18).

Modules
-------
* ``ops``      OPERATIONAL. Fermat, Miller-Rabin, Solovay-Strassen, (strong) Lucas,
               Baillie-PSW, quadratic Frobenius and AKS, each instrumented to report
               every proper gcd of an intermediate value with N.  Inputs: N, public
               parameters, a base / seed.  Never sees p, q, phi or group orders.
* ``dlog``     OPERATIONAL. Pollard kangaroo and BSGS that recover s = p+q (or q-p)
               from a Fermat (or Lucas) residue alone, then factor.  Cost ~ N^(1/4).
* ``formulas`` DIAGNOSTIC. Exact leak probabilities from the factorisation
               (order-pattern lemma and closed forms).  Uses p, q and factorisations
               of group orders.  Never imported by ``ops`` or ``dlog``.
* ``verify``   HARNESS. Exhaustive checks of the formulas and identities.
* ``experiments`` HARNESS. Leak-rate, kangaroo-scaling, AKS and partial-information runs.
"""
