(Read /Volumes/git/factoring-connection/codex/hatch_common.md and follow it.) HATCH 1, constructive + experiments.
Context: notes/crtdecode_results.md (§1.2 regime, the capacity row ~line 262: "no integer CRT capacity decoder is
known; BHKS covers polynomial ideals only"), notes/sinonobliv_results.md (~line 170: the smooth-search regime is far
below Johnson but far above capacity — lists small, no decoder), notes/pilist_results.md (H_sel, Prop JR), notes/
sihard_results.md (Thm J), subthird M2 (notes/subthird_results.md), Boneh STOC 2000 "Finding smooth integers in short
intervals using CRT decoding".
Idea to pursue: integer analogues of the polynomial tricks that reach capacity.
 (a) MULTIPLICITY: residues mod q^m are to integers what f mod (X−a)^m is to polynomials; smoothness = high valuation
     at many q, i.e. a multiplicity-weighted agreement. Formulate smooth-number search (Boneh intervals; and the 2-D
     quadratic norm version M2) as decoding a multiplicity-CRT code; compute rate/agreement in that formulation; check
     whether a Guruswami–Sahai–Sudan weighted decoder with valuation weights moves past Johnson and by how much.
 (b) FOLDING: over a Galois number field K, use the Galois action on O_K/𝔮 (Frobenius orbits, or an automorphism
     σ of K) as the analogue of folded-RS's x ↦ γx; build a "folded CRT code" over O_K whose codewords encode
     candidates and whose agreement encodes smoothness at split primes. Is there a Parvaresh–Vardy / Guruswami–Rudra
     style decoder (auxiliary relation between a and σ(a)) for it?
 Implement the most promising toy version and measure: per-witness work vs plain sieving / lattice sieving at
 40–96 bits, and the list-size / agreement curve vs Johnson. Give the exponent ledger if it scaled.
