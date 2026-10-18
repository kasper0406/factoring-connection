"""ecmavg: an explicit ECM family as a prefilter for deterministic average-case factoring.

notes/ecmavg_results.md is the write-up; notes/ecmavg_prereg.md the frozen pre-registration.
Modules:
  exponents  exact-fraction cost/exception bookkeeping (fallback bound (*), optimisation)
  curves     operational ECM prefilter over Z/N (sees only N): fixed Suyama / control
             Montgomery families, x-only ladder, staged gcds (one gcd per prime power)
  labels     DIAGNOSTIC (factor-aware): exposure labels k_p(j) computed modulo known primes,
             vectorised with numpy; never imported by operational code
  experiments, analyze   the pre-registered runs and their statistics
"""
