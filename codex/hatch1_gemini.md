(Read /Volumes/git/factoring-connection/codex/hatch_common.md and follow it.) HATCH 1, theory + literature.
Context: notes/crtdecode_results.md (capacity row ~262), notes/sinonobliv_results.md (~170–190), notes/pilist_results.md,
notes/sihard_results.md (Thm J), notes/subthird_results.md (M2), notes/multcopp_results.md ("multiplicity over Z is like
differentiating whole numbers"). Question: why do capacity-achieving list decoders (folded RS, multiplicity codes,
Parvaresh–Vardy, Guruswami–Rudra, BHKS ideal framework, Kopparty's multiplicity codes, recent plain-RS capacity
results Brakensiek et al. 2609.08005 / Jeronimo 2609.05870) exist for polynomial codes but not for integer CRT codes?
Identify the exact algebraic ingredient each uses (derivatives / Hasse derivatives, the automorphism x ↦ γx,
subspace-evasive sets, Frobenius) and whether an integer or number-field analogue exists: p-adic expansions mod q^m
as multiplicity; Galois automorphisms of O_K or Frobenius on O_K/𝔮 as folding; Fermat quotients / p-derivations
(Buium's arithmetic differential equations, δ_p(x) = (x − x^p)/p) as an arithmetic derivative! Assess whether
Buium-style p-derivations could give a CRT "multiplicity code" over Z with a decoder past Johnson in the smooth-search
regime (agreement fraction ~ 1/π(B), rate ~ F_{1−β}/(x π(B))). Give the exponent ledger if it worked, the precise
missing lemma, and a toy test (python3) if feasible. Use only read-only shell commands plus python3 in your cwd.
Final answer ≤ 2000 words, plain text.
