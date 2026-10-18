"""Agent `tower128`: what determines s = p + q mod 128 when N = pq = 1 (mod 8)?  (NEXT_STEPS item 2.)

Everything in this package is DIAGNOSTIC / HARNESS code: it knows p and q and checks mathematical claims by
exhaustive enumeration.  Nothing here is operational (rule 2 of notes/AGENT_RULES.md).

Modules
  formula     exact formula for uv mod 8 via 2-adic Dirichlet characters; exhaustive check
  classgroup  exact brute-force 2-Sylow structure of Cl(D) (primitive forms), validated against src/repair
  run_tower   computes the extended key (2-power ranks for 8 public discriminants) for all pairs, JSONL output
  analyze     ambiguity tables and conditional mutual information with the hidden bit
  symbols     quartic / octic residue symbols and Redei-symbol identities used in notes/tower_results.md section 3
"""
