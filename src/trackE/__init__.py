"""Track E: partial-information formats -> lattice (Coppersmith/Howgrave-Graham) decoders.

Module roles (AGENT_RULES §2):
  coppersmith.py  OPERATIONAL decoder: sees N + a public hint only.
  formats.py      OPERATIONAL reductions of derived formats to the affine hint.
  oracles.py      HARNESS / DIAGNOSTIC: knows the factors, produces hints.
  info_floor.py   HARNESS: exact per-modulus information-floor census.
  exp_*.py        HARNESS experiment drivers.
Operational modules never import oracles.py / instances.py.
"""
