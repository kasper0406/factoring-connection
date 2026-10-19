# rkrank_1 — blind research round on Hypothesis RK (GF(2) rank of NFS relations)

Repo (read-only for you): /Volumes/git/factoring-connection. Context: notes/hypD_results.md §5 (Design V,
Hypothesis RK, Theorem V, (B7), E3 in §7.3, codex hypD_2 review in §9.2), notes/nfsfruit_results.md (LV
structure, Lemma H, Prop R), src/hypD/rank.py (toy E3). Lee–Venkatesan = arXiv 1805.08873. Be rigorous;
label every claim proved / conditional / heuristic / false. Do not overclaim. Answer ≤ 3500 words.

## Problem
N odd composite, small factors removed (2^ω = L^{o(1)}). Field A = (m, f), m ~ N^{1/(d+1)}, f = f̂_{N,m} +
(x − m)R with R random of height H (LV-style), fixed odd d ≥ 9, factor base B = L[1/2, a]. Relations:
coprime (a, b) in a box of size L^{O(1)} with a − bm and F(a, b) both B-smooth (i.i.d. from the kept law of
one field, or the whole smooth set X_{m,f} of the box). Row of a relation: rational parities (sign,
v_p(a − bm) mod 2, p ≤ B), algebraic parities v_𝔭(a − bα) mod 2, LV §6 quadratic-character columns.
R = row space. W = {u : (u, 0, 0) ∈ R} (rational-only part; genuine-certificate version W_true after LV §6).
Λ/2 = {u : ∏ p^{u_p} ∈ (Z/N)^{×2}}. Known: W_true ⊆ Λ/2. Vallée vectors (parity vectors of smooth
x² mod N, |x² mod N| ≤ N^{2/3}) lie in Λ/2.

**Hypothesis RK.** With probability ≥ L^{−o(1)} over the field, codim of W_true ∩ U in U is ≤ L^{o(1)}, where U
is (Λ/2 restricted to well-covered primes, or the span of the Vallée target supports). Then L^{o(1)} Vallée
relations force a dependency containing a Vallée row, and Dixon/Vallée fibre symmetry makes it fruitful,
giving rigorous L[1/2, √(2/3)] (Theorem V, conditional also on LV relation supply at L[1/2] and LV §6).

## Tasks (attempt a proof first; if you cannot, find the sharpest precise obstruction)
1. **Dual obstructions.** W^⊥ = proj_rat(R^⊥). A dual obstruction is (w, y) with w(a − bm) = y(a − bα) for
   every relation. Classify the forced ones: the ω characters (χ_p, χ_𝔓) for p | N; archimedean ones (sign
   of a − bm vs sign at the real root near m); leading-coefficient / |S|-parity columns; uncovered or rarely
   covered primes; anything else (Schirokauer-type maps, small-conductor characters, units, Cl(K)[2]).
   Is "W^⊥ = Gen + archimedean + sparse-on-rare-primes" plausibly true? Provably?
2. **Proof routes.** Try: (i) uniform rank growth — for every subspace V of non-forced obstructions,
   Pr_relation[relation ∉ V^⊥] ≥ dim V/(D·L^{o(1)}) (D = #columns), which gives codim ≤ r after D·L^{o(1)}
   relations by a coupon-collector sum; (ii) Fourier/Parseval ("large sieve over GF(2) functionals")
   counting of near-obstructions via collision probabilities of sums of relations; (iii) Halász-type mean
   values of ±1 multiplicative functions on smooth numbers (de la Bretèche–Tenenbaum friable
   Turán–Kubilius etc.) — what do they give on the rational side alone, and why does the joint
   rational/algebraic constraint break them; (iv) Hafner–McCurley / Buchmann random-exponent relation
   lattice generation: can any cheap NFS-sized relation carry random exponents (exact or quasi-uniform
   on Λ/2 given the transcript)? (v) many independent fields (W = Σ_i W_i).
   For each: where exactly does it fail, and what is the minimal extra hypothesis?
3. **Barriers.** Is there a rigorous reason why averaging over the field (≤ N^{O(1)} fields) cannot give
   RK via union bounds over 2^D functionals (a "fixed-test" barrier)? Is RK plausibly FALSE anywhere
   (e.g. structured fake obstructions of large dimension)? Give explicit families if you find any.
4. **Relaxations.** Show (or refute) that a weaker codimension bound D^κ suffices: with a sampler of exponent
   θ, Design V costs L[1/2, √(2θ/(2 − κ))] if codim ≤ D^κ = L[1/2, κa]. Check this formula.
5. **Literature.** Stange 2022 (arXiv 2211.06821) Hypothesis 3.1 and what she proves; Hafner–McCurley 1989,
   Buchmann 1990, Schnorr–Seysen–Lenstra / Lenstra–Pomerance 1992 (which parts are rigorous, which rank
   arguments they use); Balister–Bollobás–Morris 2018; rank of sparse random GF(2) matrices (Costello–Vu,
   Coja-Oghlan et al.); any rigorous span/generation statement for sieve relation matrices (2000–2026).
   Has anyone rigorously beaten Lenstra–Pomerance c = 1?
Give a short final verdict: what is proved, the best precise reduction of RK, and the single most promising
next step.
