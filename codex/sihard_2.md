# sihard round 2 (ADVERSARIAL REFEREE): hardness of constructing smooth numbers in short intervals near kN

You are an adversarial referee (analytic number theory of smooth numbers, Coppersmith / Howgrave-Graham / Boneh /
Guruswami–Sahai–Sudan CRT list decoding, lattice algorithms, index calculus for discrete logarithms, complexity theory:
NP-hardness, TFNP, ETH). Read the repository file `notes/sihard_results.md` (agent sihard) in full, and the code
`src/sihard/core.py`, `src/sihard/reduction.py`, `src/sihard/consequences.py`, `tests/test_sihard.py` (you may run
`.venv/bin/python -m pytest -q tests/test_sihard.py` from the repository root). You may consult
`notes/suffmap_results.md` §1.4 (Prop SI), `notes/pilist_results.md` §2.2 (Prop JR), `notes/siconstruct_results.md`
(a parallel agent) and the primary literature via web search. Work in the scratch directory for computations.

Your job is to BREAK the note: find wrong statements, mislabelled statuses (PROVED / model-PROVED / KNOWN /
HEURISTIC), quantifier errors, hidden costs, missing hypotheses, wrong citations, and models that are too weak to
mean anything. Focus in particular on:

1. **Theorem J** (§5.2). Check every step of the unweighted proof: the averaging step, the divisor decomposition,
   the uniformity of Hildebrand's Ψ(Y, y) ≤ C₀Yρ(u) bound over 1 ≤ Y ≤ 3X at y = L[1/3], the Nicolas–Robin divisor
   bound, the injection/Rankin bound (b), Case A's concavity argument, Case B's π^{0.7} = o(σ) estimate, the
   small-region case, and the final o(1). Then the weighted sketch: the Cauchy–Schwarz class reduction, tail
   truncation, the good/bad class bookkeeping (threshold R₀ = C n^{4/7} log n; is the bad r-mass really o(t)?), and
   the per-class inequality. Then extensions (a) averaging over N, (b) Herrmann–May batched regions, (c) families of
   moduli. Is the class 𝒥 correctly specified (does Boneh's Theorem 3.1 really fit; GSS's actual theorem)? Is the
   "sharpness" claim right? Is the reading (Johnson class ≡ divisor trick at leading order; escape needs
   ≈ (2/(3κ))F_{2/3} bits of x-dependence) justified? Is the result already known (Boneh's §3.1 counting bounds,
   Coppersmith–Howgrave-Graham–Nagaraj on divisors in residue classes, Cohn–Heninger, anything)?
2. **Prop DLI** (§2.2): the descent with modulus M = ⌊jp/q⌋, the size of v, whether L_M = L_p to leading order is
   enough, the ECM booting constant 3^{1/3} (is it the right classical constant; is the √(4c/3) ECM cost right; are
   there better boostings — early abort, Barbulescu, Guillevic, Joux–Lercier — and what constant do they give?),
   whether the constructor could also boot, the linear algebra mod p − 1, and whether the conclusion "≤ L_p[1/3, 1.442]"
   is correctly labelled. Is anything like this in the literature?
3. **Lemma BB, Prop S, Prop MM, the Wagner bit-budget argument, the lattice dimension lemma** (§§3, 5.1, 5.3–5.5):
   are the models meaningful, are the proofs right, and does the collapse table (§5.6) overstate anything?
4. **§6**: Prop NC (NP ∩ coNP for h ≤ poly(log x): is the coNP certificate really polynomial; what about y-smoothness
   testing being factoring-hard?), Prop TF (Megiddo–Papadimitriou applies to Cook reductions? to promise problems?),
   the list of all-interval existence results (Balog 1987, Harman 1991, Matomäki–Radziwiłł 2016 Cor. 1,
   Friedlander–Lagarias 1987 — check the exact statements), the ETH blow-up claim, Lemma FP, and Prop XS (check the
   reduction's inequalities in `src/sihard/reduction.py`, the parsimony, the NP-membership with unbounded exponents,
   and the #P remark). Is there a known NP-hardness result for smooth numbers in intervals that we missed?
5. **§7 barrier statement and §0 verdict**: overclaims? Is "no barrier is known" accurate? Anything missing (e.g.
   other natural models: adaptive lattice/decoder choice, algebraic constructions, polynomial identities à la
   Friedlander–Lagarias, number-field norms, smooth values of polynomials evaluated at points near (kN)^{1/d})?
6. **§8 experiments**: are the post-hoc E1 diagnoses sound? Is there a known explanation of a modulus-dependent bias
   of smooth counts along {kN} at toy sizes?

Output: numbered findings with severity FATAL / MAJOR / MINOR / NIT, exact location, the problem, a concrete fix. End
with an overall verdict (ACCEPT / ACCEPT WITH FIXES / MAJOR REVISION / REJECT).
