"""avgtime: complete deterministic factoring with average-cost guarantees on S_n.

notes/avgtime_results.md is the write-up.  Modules:
  exponents  exact-fraction audit of Theorem D' constraints and exponent optimisation
  detector   operational small-order detector (BSGS for min ord_r(2) <= T, r | N)
  ladder     operational ladder: detector, DCL-det at growing H, complete fallback
"""
