# equivtarget: treat "secretly ≡ FACT / ≡ DL" as TARGETS, not kills (common preamble)

You are a research consultant on /Volumes/git/factoring-connection (do NOT edit it; scratch only in your cwd).
User's only goal: factoring general N in L_N[α], α < 1/3 (heuristic fine; constants worthless).

## The reframing (the user's point)
Throughout phases 18–32 many routes were dropped because the needed primitive is "factoring-equivalent",
"FACT ∧ DL", "would solve prime-field DL at L[β]", or "poly-time factoring in disguise" (e.g. suffmap Thm SS,
navigate Prop H / Thm DLB, prepread Cor O1-N, sihard, suffmap Prop DLT, siconstruct Prop CV (promise-CVP),
smoothbuild, ecount/tracedata (#E(Z/N)), genusleak (class groups), factorbits (prime-base addition ≡ FACT; Zech ≡ DL),
subthird M2 (class-group DL), selfred, leakmap). That is a kill only if factoring/DL are hard — exactly what we are
betting against. A reduction "T ⇒ FACT" says nothing about T's native difficulty; T may have structure or
special easy regimes that factoring lacks. Phase 29 (notes/native_results.md, nativeinv/nativemod/nativeana) already
did this for 55 *invariants* (φ, σ, class numbers, a_f(N), …): none beats factoring natively; a single a_f(N) is open.
It did NOT revisit the *constructive primitives and DL-type targets* above. Do not redo phase 29.

## Key question
For each such target T: is there a regime, special family, or transfer where T is KNOWN to be easy (or much easier
than generic), and can the instances that factoring produces be steered into that regime cheaply? Examples of easy
regimes in the DL / isogeny world (check whether any analogue applies):
- small-characteristic finite-field DL: quasi-polynomial (BGJT, Kleinjung–Wesolowski), via Frobenius structure;
- medium-prime / composite-degree fields: L[1/4]-type or better tricks (Joux), tower/exTNFS;
- anomalous elliptic curves (Smart–SSSA, p-adic log): polynomial; MOV/Frey–Rück pairing transfer to F_{p^k};
- Weil descent / GHS; index calculus on high-genus curves (L[1/2]→ better in some families);
- SIDH broken in poly time via Kani's lemma / higher-dimensional isogenies given auxiliary torsion images;
- hidden-structure attacks: lattice/Coppersmith with extra hints; HNP; LLL on planted structure.
Also: DL-equivalence to F_p-DL at L[β<1/3] is itself an open target — if a reduction routes FACT to DL in a group
with an easy regime, that is the jackpot.

Output (≤ 1800 words, plain language; label proved/KNOWN/heuristic/speculative; cite file:line):
- candidates ranked by promise, each with: target T and the reduction (direction, cost blow-up, does T at L[α]
  give FACT at L[α]?), the easy regime/transfer, the steering step needed and its cost, the precise missing lemma,
  a cheap toy falsification test (pre-register before running; run if feasible), P(α<1/3 heuristic).
- explicit kills of your own ideas.
