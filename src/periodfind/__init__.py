"""periodfind: classical non-black-box period finding for x -> a^x mod N (phase 17).

Modules:
  observables -- OPERATIONAL observables g(y, N) (functions of the residue y and N only)
  spectra     -- HARNESS: orbit enumeration, Fourier spectra on Z/r, metrics, transfer-identity check
  sft         -- OPERATIONAL sparse Fourier (KM low-bit bucketing on Z_{2^m}) + continued fractions
  generic     -- HARNESS: generic-ring decision trees, p/q runs, split rate, visible mass
  run_experiments -- HARNESS: E1-E4 of notes/periodfind_prereg.md
"""
