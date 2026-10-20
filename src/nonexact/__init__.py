"""nonexact (phase 17): O1, non-exact / evaluated compressions of split polynomials modulo N.

Modules
- core.py: exact helpers (polynomials, reductions mod (N, F), lattices via flint LLL, exact
  Fincke-Pohst enumeration, the linearised-lattice counts of Prop LIN, Prop IH relation extraction).
- run_experiment.py: the pre-registered experiments E1-E3 (notes/nonexact_prereg.md).

All code here is DIAGNOSTIC/HARNESS unless a function is marked OPERATIONAL. Operational functions
receive only N, public parameters and a seed (rule 2).
"""
