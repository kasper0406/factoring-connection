"""pm1bsgs: evaluation of the user's "p-1 stage 1 + BSGS on b^N = b^x (mod p)" idea.

Modules
-------
operational : OPERATIONAL algorithms (see only N, public parameters, seed):
              stage 1, BSGS stage 2 (mode "N" = user's hybrid, mode "order"),
              prime-only stage 2, birthday stage 2, plain BSGS.  Instrumented.
costmodel   : OPERATIONAL-SAFE exact cost formulas mirroring operational.py
              (parameters only, no factor knowledge).
harness     : DIAGNOSTIC / HARNESS (knows p, q): smooth parts, element orders,
              exact success predictors, exact failure probabilities.
              Never imported by operational.py or costmodel.py.
experiments : experiment drivers (harness side).
"""
