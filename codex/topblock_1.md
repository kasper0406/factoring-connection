# topblock_1: blind round. DOES THE FLOOR-SUM LADDER'S TOP-BLOCK WINDOW ORACLE GIVE ANYTHING BEYOND KNOWN-TOP-BITS FACTORING?

You are a research mathematician (computational number theory, lattice reduction / Coppersmith
methods, cryptanalysis) with web search. The repository at /Volumes/git/factoring-connection is
read-only context (an AI-assisted factoring research project). Write scratch files only under your
current directory. You may read notes/analytic_results.md §3 (Theorem L_j), notes/floorsum_results.md
§6 and notes/floorhard_results.md §4.2, and cite them instead of re-deriving.

Only CLASSICAL algorithms matter. Be honest: a precise obstruction is as valuable as an algorithm.
Label every claim PROVED (give the proof), KNOWN (give an exact reference you have checked: authors,
venue, year, theorem number), HEURISTIC, CONJECTURED, or UNVERIFIED.

## Setting

N = pq balanced (p, q ≈ √N). Fix j ≥ 2, 0 < c ≤ ½, B = ⌈N^{1/(j+2)}⌉. A block is [x, x+H) with
H = ⌊cx/B⌋; at the top (x ≈ √N) H ≈ c·N^{1/2 − 1/(j+2)}. P_j(h) = Σ_{i=0}^{j} (−1)^i N h^i / x^{i+1}
(Taylor polynomial of N/(x+h), denominator D = x^{j+1}); remainder |R(h)| ≤ ρ = N H^{j+1}/x^{j+2} < 1.
The ladder's oracle answers, for any rational shift and any subinterval I of [0, H), the count
#{h ∈ I : {P_j(h)} ∈ window} (equivalently, exact floor sums Σ⌊P_j(h) + s⌋). Divisors x+h | N land in
the window [0, ρ] (even j). With a polylog oracle the ladder factors balanced N in about N^{1/(j+2)}
oracle calls. Known: (a) the oracle for ALL x is NP/#P-hard, but only via blocks at x ≤ N^{3/8};
(b) at top blocks, exact-value queries are Coppersmith-easy; (c) the best exact floor-sum oracle
costs H·(N/x³)^{1/3} per block, giving N^{1/3} overall; (d) a bivariate Coppersmith lattice for the
window congruence r ≡ num(h) (mod x^{j+1}), 0 ≤ r ≤ ρx^{j+1}, resolves only sub-blocks of length
≈ x·N^{−1/3}. Locating p in a top block is "factor N given the top (1/(j+2))·log₂N bits of p".
j = 2 is Coppersmith's half-the-bits theorem (N^{1/4}); j ≥ 3 would beat it.

## Tasks (answer each; be precise)

1. **Equivalence.** Is "top-block window oracle at degree j" exactly equivalent to "factor N from
   (1/(j+2))·log₂N top bits of p", or does the window structure carry extra information
   (correlations across consecutive blocks, the Taylor structure, multiple windows/shifts, windows of
   several degrees j)? Give reductions in both directions with their exact costs and any loss
   (e.g. near-divisors / candidate counts). If there is extra information, quantify it.
2. **Partial-key-exposure frontier.** State precisely (with references you checked): Coppersmith's
   high-bits result; the cost with fewer bits (brute force over missing bits); May's survey
   "Using LLL-reduction for solving RSA and factorization problems" (2010) Theorem 1; Herrmann–May
   (non-consecutive known bits / multiple unknown chunks) and its bound as the number of chunks grows;
   May–Ritzenhofen implicit factoring; Boneh–Durfee–Howgrave-Graham for p^r q; Heninger–Shacham.
   Is factoring from fewer than half of p's top bits in polynomial time known to be impossible for
   lattice (Coppersmith-type) methods? State any known lower-bound / optimality arguments precisely
   (determinant/dimension arguments, root-count/list-size arguments, "optimality of lattices for the
   Coppersmith technique", the Johnson-bound analogy for CRT codes) and their exact scope.
3. **Multi-block / multi-window lattices.** The ladder supplies one window condition per block, all
   about the same p, plus many shifts and subintervals per block. Can any lattice combining several
   blocks' or several windows' conditions (mod the known D = x^{j+1}, mod the unknown p, or over Z)
   beat the single-block Coppersmith range N^{1/4}, in the way Herrmann–May combine several chunks of
   known bits? Cost it (asymptotic determinant/dimension analysis) and say what toy experiment would
   decide it.
4. **Verdict.** Either a mechanism beating N^{1/4} (or N^{1/3}) via structured windows, or a precise
   statement that the ladder is exactly as hard as partial key exposure with a given fraction of bits,
   with where that barrier is believed to sit.

Keep the answer under ~4000 words. Put the verdict first.
