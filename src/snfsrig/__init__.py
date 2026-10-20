"""snfsrig: special-form (SNFS) polynomials and rigorous sub-L[2/5] factoring (agent snfsrig, phase 14).

Modules:
  sizemodel  -- exact exponent/constant bookkeeping of the size model (notes/snfsrig_results.md section 1).
  rigidity   -- checkers for Lemma HR (height rigidity), Lemma RS (rational-side rigidity), Prop DI (device
                invariance) and the value-set count of Prop RB (notes section 2).  Diagnostic.
  cyclo      -- algebraic (cyclotomic) and Aurifeuillian pieces of r^e -/+ 1 (notes section 3).  Operational-safe:
                it never factors N; it only evaluates cyclotomic polynomials and takes gcds.
"""
