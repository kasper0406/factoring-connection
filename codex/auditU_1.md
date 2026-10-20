You are a blind, adversarial referee (analytic number theory + computational number theory). The repository at
/Volumes/git/factoring-connection is read-only for you; write any scratch scripts under your current directory.

## The claim to referee

`notes/auditU_blind_claims.md` (a verbatim excerpt of an internal note, with review sections removed) claims
**Theorem U**: a rigorous Las Vegas algorithm that completely factors every integer N in expected time
L_N[2/5, c] with a finite constant c. This would be a major result: the best rigorous bound we know of for arbitrary N
is Lenstra–Pomerance 1992, L[1/2, 1+o(1)]; heuristic NFS is L[1/3]. The claim is conditional on an internal,
unrefereed fixed-degree chain whose cleanest form is the preprint `paper2/main.tex` (cited as [pp]), plus
`notes/lem85_results.md` (Lemma 8.5′, a restricted Bombieri–Vinogradov for smooth numbers) and
`notes/recurse_results.md` (Lemma S″, Lemma F, Lemma OS, the recursive shape test; §3 is the size model).
Theorem U re-runs that chain at **growing** degree d ≍ D^{1/5} and k ≍ D^{2/5} (D = log N/log log N), with B = L[2/5],
special-q primes L[3/5], and claims every precision loss is polylogarithmic ("τ = 0"), so that the split success
probability (γ/(16 log N))^{2k} is L[2/5]^{O(1)}.

**Blindness rules.** Do NOT read STATUS.md, notes/lv5unif_*, codex/lv5unif_*, notes/*audit*.md (except
notes/auditU_blind_claims.md), codex/auditRA_*, codex/thmA_*, or any other file that records a review verdict on this
claim. You may read paper2/main.tex, notes/lem85_results.md, notes/recurse_results.md, notes/{pcest,ldbound,lv6,hypD,lemmaM}_results.md,
src/lv5unif/*.py, and the web (Lee–Venkatesan arXiv 1805.08873 §§5–6, 8; Hildebrand–Tenenbaum 1986; Drappeau 2015;
Harper 1208.5992; Stark 1974; Lagarias–Odlyzko 1977; Gordon 1993; anything else).

## What to do

Try hard to break Theorem U. Assume the fixed-degree chain [pp] is correct *at fixed d* (it has been audited at fixed d), but
check every place where Theorem U says "[pp]'s proof, with these explicit inputs" or "uniformity checked": a hidden
constant that is exponential in d (2^d, d^d, C_d) or a (log N)^{f(d,k)} with f unbounded anywhere in the relation-supply
probability γ is fatal, because γ is raised to the power 2k ≍ D^{2/5}. In particular:

1. **The four τ = 0 repairs.** (a) b-window ((1−1/d)V, V] with one radius V for all special-q primes; (b) envelope
   X_F = (2d+1)MHV^d; (c) Lemma Kern-U (triangular fibre count, Pr ≤ e/((2H+1)²(bb′)^{d−2}) for bb′ ≤ (2H+1)/d, used with
   H ≥ dV²); (d) Lemma 8.5′ in the level-gap regime (R′) with Y = L[3/5]. Is each repaired estimate true, with only
   poly(d) loss? Does the narrower window or the common V break the first moment (FM-U), the bad-modulus bound, the
   special-q lattice search ([pp] lem:reg), or Theorem K-U? Does H ≥ dV² change θ?
2. **Uniformity at growing (d, k).** Go through every N ≥ N₀(·), every C_d and every o(1) in FM, Theorem R, Lemma W, L2_fam,
   LD*, Theorem 4.3/D, Lemma S″, Theorem K, LV6 and Lemma 8.5′ in [pp]/lem85, and decide whether each holds with d, k growing.
   Specific checks: Hildebrand–Tenenbaum / de Bruijn / CEP error terms at u ≍ D^{2/5}; Drappeau Lemme 5, Harper, the large
   sieve at y = B = L[2/5] and moduli ≈ MV; Tenenbaum's Ψ_q theorem; the computable S (Lemma CS: is m = S²ρ_rρ_F really
   polylogarithmic?); η₀ = D^{−1/2}; Lemma IR (irreducibility) and all other field rejections at growing d; whether any proof in
   [pp] silently uses d fixed (e.g. compactness, "O_d(1)", a union bound over 2^d objects, or a fixed number of base digits).
3. **LV6 at growing d.** Character primes of (2d)!·poly(log N) bits via Stark's Lemma 8 + Lagarias–Odlyzko; the exceptional-class
   lemma; square roots in degree-d fields of exponent vectors of size L[2/5]; costs.
4. **Recursion at Λ_{2/5}.** Lemma F_α; the recursive shape test on |u₀|, |v₀| ≤ √(2N); oracle substitution.
5. **Cost.** Re-derive c = max(2a, a+s+max(a,σ), (a+s)/(1−2^{−2/5})), s = 2κ(1+C_γ), and the geometric argument (§13)
   for α = 2/5. Is c finite without knowing C_γ, i.e. is C_γ really independent of N, d, k?
6. **Sanity.** If this works, why did nobody do rigorous L[2/5] before (Lee–Venkatesan already prove rigorous L[1/3]
   congruences of squares)? Is there any known obstruction? Compare with Gordon 1993's discrete-log descent and with LV §8.

## Output

Write your report in markdown (it is captured automatically). Start with a one-line whole-theorem verdict:
ACCEPT / ACCEPT WITH FIXES / MAJOR GAP / FALSE. Then a table of issues (FATAL / MAJOR / MINOR, location, why, fix if any),
then the single **weakest link**, then answers to items 1–6. Be concrete: quote the line, give the counterexample or the
exact estimate that fails. Do not pad; do not accept "line by line as [pp]" without checking the [pp] proof yourself.
