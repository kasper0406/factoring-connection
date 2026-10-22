# sinonobliv round 2 (adversarial review)

You are an adversarial referee (analytic number theory, lattice decoding, factoring). The repo
/Volumes/git/factoring-connection is read-only for you; scratch goes in your current directory; python with sympy,
gmpy2, python-flint, numpy is at `/Volumes/git/factoring-connection/.venv/bin/python`. Web search is available.

Review `notes/sinonobliv_results.md` (version 1, commit 04a1f52) together with its code `src/sinonobliv/*.py`,
tests `tests/test_sinonobliv.py`, pre-registration `notes/sinonobliv_prereg.md` and raw data
`experiments/raw/sinonobliv_*.jsonl`. Background (take as given, but check how the note uses them):
`notes/sihard_results.md` §5.2 (Theorem J and its proof), `notes/siconstruct_results.md` §2.2 (Lemma SR/RC),
`notes/sicheck_results.md` §2.1 (Lemma I, Lemma CP), `notes/suffmap_results.md` §1.4 (Prop SI). Your own blind round is
`codex/sinonobliv_1.out`.

Try hard to break things. In particular:
1. **Prop JC** (Theorem J for fixed chains): is the reduction to Theorem J at scale (X/D, V/D) valid (uniformity of
   Theorem J's o(1) at X′ = X^{1−o(1)} with the same y; the D > V case; weighted case)?
2. **Lemma IC / Prop RD**: is "no informational barrier can exceed L[1/3, κ]" correctly stated? Is the CRT-code
   reformulation (weights, prime powers, Johnson radius, list size) right? Is the decoder threshold a < 1.1416 right
   (c = (3/2)(4/9)^{1/3}(1 + a)^{2/3} for a ≥ 1/3)? Any literature on CRT decoding beyond Johnson / capacity-achieving
   CRT-type codes that bears on this?
3. **Prop FF** (fixed families must be dense in the smooth numbers): proof, multiplicity handling, and whether the
   "reading" (correlated families fail; dense families are SI again) overreaches. Is "only L[1/3, C] multipliers ⇒
   p ≥ L[1/3, −C] per window" justified?
4. **Prop MS** (independent windows, partial-ACD region, any y-smooth P, 2 ≤ m ≤ ln h): check every inequality of
   the two cases (j₁ = mt/(m−1) + 1; π ≤ 60t in case 2; ρ-ratio bound; Hildebrand uniformity when d is close to X).
   Is the region granted correctly (σ_m = t^{m/(m+1)}π^{1/(m+1)}), and does "certified d ≥ h + 1" follow?
5. **Prop MS′** (same N, d ≤ Kh): the conditioning on N mod d and the arithmetic-progression bound.
6. **Prop NV, Prop LC** (route 4): is Prop LC's ideal-counting heuristic sound (unit index, class number, the claim
   that the constructed side is ≍ N/(h c_K))? Is anything stronger provable?
7. **Experiments**: re-score the registered predictions from the raw data; check the post-hoc local-factor model
   (posthoc.py) and the E2B seed-spread diagnosis; check E3's certificate threshold (uncert.py) and classification;
   any code bug that changes a number.
8. **Verdict and labels**: overclaims, missing caveats, wrong attributions, places where the note says PROVED but
   only a model statement holds.

Output: a table of findings (#, severity FATAL/MAJOR/MINOR/NIT, location, finding, suggested fix), then details with
explicit counterexamples or computations where you disagree, then an overall recommendation (accept / minor / major
revision). Be concrete and terse.
