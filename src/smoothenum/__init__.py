"""Agent smoothenum (phase 15): output-sensitive enumeration of doubly-smooth NFS pairs.

target.py  : L[1/3] cost programs with a generic enumeration-cost term (analysis; floats only in minimisers)
barrier.py : the lattice-query / generation-model yield bounds (Dickman rho, analysis only)
nfsdata.py : polynomial pairs (seeded base-m, or CADO .poly files), factor bases, exact norms
mech.py    : mechanisms (box line sieve, special-q, special-(q1,q2), smooth-modulus lattices, bottom-up)
run_proto.py : pre-registered prototype experiment driver (streams JSONL, chunked, low memory)
"""
