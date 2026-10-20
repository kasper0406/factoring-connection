"""nonsmooth (phase 16, route 3): factoring paradigms not based on smooth numbers.

Module roles (AGENT_RULES §2):
  equiv.py   Equivalence map. ``reduce_*`` functions are OPERATIONAL reductions: they see N and
             the value of one invariant (as an oracle answer) and nothing else. ``oracle_*``
             functions are DIAGNOSTIC: they compute the invariant from the factorisation, to
             feed the reductions in tests. Operational code never imports an oracle.
  forms.py   Binary quadratic forms (reduction, composition, class numbers by enumeration).
             Class-number enumeration is O(D) and DIAGNOSTIC-cost only.
  cheng.py   OPERATIONAL: the CM / anomalous-curve special-purpose factorer (Cheng 2002;
             Shirase 2017) for class-number-one discriminants. Sees N only.
  census.py  DIAGNOSTIC: squarefree part of 4p-1 over random primes (coverage of the
             anomalous-CM paradigm). Uses factorisation of 4p-1 for primes p it generates.
"""
