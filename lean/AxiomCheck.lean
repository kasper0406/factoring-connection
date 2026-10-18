/- Run with `lake env lean AxiomCheck.lean`: prints the axioms each main result depends on.
   Expected: only `propext`, `Classical.choice`, `Quot.sound` (no `sorryAx`, no custom axiom). -/
import FactoringLean
open FactoringLean
#print axioms exponent_identity
#print axioms pow_exponent_identity_zmod
#print axioms pow_exponent_identity_modEq
#print axioms coverage_criterion
#print axioms coverage_criterion'
#print axioms lemma1_class_lock
#print axioms lemma2_lower
#print axioms lemma2_upper
#print axioms lemma2_index
#print axioms lemma2_window
#print axioms lemma3
#print axioms lemma4
#print axioms lemma4_exponent_difference
#print axioms lemma5_one_sided
#print axioms lemma6_clean_band
#print axioms lemma6_clean_band_prime
#print axioms coprime_dvd_mul_iff
#print axioms lemma7_core
#print axioms lemma7_full_iff
#print axioms lemma7_coprime_zmod
#print axioms clean_coprime
#print axioms lemma8
#print axioms lemma8_sharp
#print axioms card_le_of_pairwise_modEq
#print axioms lemma10_single
#print axioms lemma10_union
#print axioms lemma10Modulus_ge
#print axioms weights_eq_of_same_class
#print axioms lemma9_of_largeSieve
#print axioms exceptional_Eb_of_largeSieve
#print axioms single_target_extraction
#print axioms single_target_from_lock
