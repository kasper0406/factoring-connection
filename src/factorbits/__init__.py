"""factorbits (phase 21): prime-exponent addition vs factoring, partial information about factors.

Modules
-------
equiv    OPERATIONAL reductions that use an oracle passed in by the caller (addition / successor /
         partial-output oracles for the prime-exponent representation; Zech vs discrete log).
         They never factor anything themselves except by trial division up to a stated bound.
harness  HARNESS oracles (these factor with sympy or brute force; never imported by equiv).
pure     The pure generic-ring decision-tree model: class partition of primes, exact split
         probabilities and the bounds of Theorem PT (diagnostic; uses the factors by design).
gl       Goldreich-Levin list decoding of p from a planted <p, r> predictor (toy demonstration).
"""
