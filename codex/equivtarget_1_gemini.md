(Read /Volumes/git/factoring-connection/codex/equivtarget_common.md first and follow it.)

## Your focus (G): transfers from known easy regimes (fresh eyes)
Start from the easy-regime list in the preamble and work backwards: for each easy regime R (small-char quasi-poly DL,
Joux medium-prime / composite-degree tricks, anomalous-curve p-adic log, MOV/pairing transfer, Weil descent/GHS,
Kani / higher-dimensional isogeny attacks given torsion images (SIDH break), Coppersmith/HNP with hints, others you
know), ask: is there a factoring-equivalent (or FACT∧DL) problem that, for the N we are given, can be cast into regime
R at subexponential cost below L[1/3]? Think specifically about:
- formulations of factoring as DL / order / isogeny / torsion problems (e.g. order of (Z/N)^× elements, #E(Z/N),
  Cl(−4N), isogenies between curves over Z/N, Jacobians of curves over Z/N, Weil restriction of Z/N-objects);
- what "auxiliary torsion information" or "Frobenius structure" the composite modulus N provides for free, and
  whether a Kani-type gluing or a BGJT-type Frobenius trick could exploit it;
- why previous repo kills (factorbits PT, FROB/PRIM barriers, navigate DLB, selfred SR1/SR2) do or do not apply to
  your transfer — read those notes.
Give ≤ 3 surviving candidates with ledgers, missing lemma, toy test (run with python3 if cheap), P(α<1/3).
Use only read-only shell commands plus python3 in your cwd. Print the full answer as your final response.
