"""ptebarrier: barrier side of route 1's PTE escape (systeq section 6).

DIAGNOSTIC / HARNESS code only: exact integer arithmetic, no factoring of N anywhere.

- rigidity.py : Lemma R thresholds and a step-by-step certificate of its proof (archimedean
                rigidity of split polynomials whose (multiplier-twisted) difference is small and
                of low degree), the View-A support lemma (Prop S), and the "one presentation,
                one relation" invariant (Theorem W).
- enum.py     : exhaustive enumeration of split polynomials with roots in [-B, B].
- run_experiment.py : pre-registered runs E1, E1b, E2, E3 (notes/ptebarrier_prereg.md).
"""
