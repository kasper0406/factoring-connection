# rkrank_2 — research round on rigorous samplers of small modular squares (θ < 2/3)

Repo (read-only for you): /Volumes/git/factoring-connection. Context: notes/hypD_results.md §5 (Design V,
Theorem V) and notes/nfsfruit_results.md §5. Be rigorous; label claims proved / heuristic / false. ≤ 3500 words.

## Setting
A *fibre-symmetric small-square sampler* of exponent θ is a randomized algorithm that, on input N (odd
composite, no prime factor ≤ L^{o(1)}), runs in time L^{o(1)} and outputs x with |x² mod N| ≤ N^{θ+o(1)} such
that for every output x and every ε ∈ μ₂(Z/N), Pr[εx] ≥ Pr[x]/κ with κ = L^{o(1)} (e.g. quasi-uniform on an
ε-closed set such as {x : |x² mod N| ≤ h}), and the value z = x² mod N is B-smooth with probability
≥ ρ(size)·L^{−o(1)} (e.g. via quasi-uniformity of z over [−h, h] ∩ square classes plus Pomerance-type
counting). Vallée (Math. Comp. 56 (1991) 823–849) gives θ = 2/3.

## Tasks
1. **Check the consequences.** (a) Dixon-style with such a sampler: c = √(2θ) (so θ < 1/2 alone would beat
   Lenstra–Pomerance c = 1 with no rank hypothesis). (b) Design V (cheap LV rows + L^{o(1)} sampler rows,
   with RK): c = √θ, needing fixed d with 3/(d+1) ≤ θ/2. Verify both, including hidden losses (2^ω, κ,
   smoothness counting in square classes for z ≤ N^θ via products of two smooth numbers, Pollard–Strassen
   preprocessing), and whether the cost per sample may be L[1/2, s] for small s > 0 instead of L^{o(1)}.
2. **Why 2/3 in Vallée.** Reconstruct the mechanism (cells x₁ + t, |t| ≤ √h, linearised 2D lattice
   s ≡ x₁² + 2x₁t (mod N), area h^{3/2} ≥ N; Farey handling of degenerate cells). State precisely where 2/3
   comes from.
3. **Attempts below 2/3.** Try, and for each give either a proof or the precise failure:
   (a) higher-dimensional Vallée (x = x₁ + Σ t_i w_i; I believe smallness of the quadratic part forces a
   rank-1 collapse — prove or refute); (b) Coppersmith / Howgrave-Graham / Jochemsz–May for small values of
   (x₁ + t)² mod N with |t| > √h (compute the lattice condition; compare with the needed T·h ≥ N);
   (c) Elkies-style points near curves (ANTS 2000) for x² − kN small; (d) ratio forms x² ≡ u/v with
   |u|, |v| ≤ N^{θ/2}; (e) class groups / binary quadratic forms (Lenstra–Pomerance, Schnorr–Seysen, Hafner–
   McCurley random-exponent relations: they behave like θ = 1/2 objects — make this precise: is there a
   "relation sampler of exponent θ" interpretation of LP, and could a θ < 1/2 analogue exist there?);
   (f) higher-degree "small squares": x with x² ≡ g(t) for a polynomial with small values; (g) anything
   else you know (Schnorr's CVP factoring claims, Ajtai-type sieving, Hittmeir, Harvey, Boneh's CRT
   decoding, Coron, Cohn–Heninger).
4. **Barriers.** Can you prove that any sampler of a natural class (cells = 1-parameter families solved by
   2-dimensional lattice enumeration; or by Coppersmith with polynomially many shifts) has θ ≥ 2/3? Is there
   a reduction showing that a fibre-symmetric sampler with θ < 1/2 (or θ < some bound) would itself imply
   fast factoring or contradict something known?
5. **Literature 1991–2026.** Search for any improvement of Vallée's 2/3, any rigorous quasi-uniform
   small-square or small-residue generator, and any rigorous factoring with c < 1. Give exact citations.
Final: short verdict (is any θ < 2/3 known/provable; best barrier).
