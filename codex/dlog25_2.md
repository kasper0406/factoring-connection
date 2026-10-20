You are an adversarial referee (computational + analytic number theory). The repository at
/Volumes/git/factoring-connection is read-only for you; write scratch scripts under your current directory.

## The claim

`notes/dlog25_results.md` claims **Theorem DL25**: a Las Vegas algorithm computing discrete logarithms in F_p^* for every
prime p in expected time L_p[2/5, c] with c finite, "proved from U", i.e. conditional only on the internal unrefereed chain
behind the project's factoring Theorem U (`notes/lv5unif_results.md`, audited in `notes/auditU_audit.md`) and on Theorem U
itself, used as the shape-test subroutine. Its DL-specific parts: Lemma S_p″ (split counting for T uniform on F_p^*),
Lemma T_p-U (hidden labels with a randomized factoring subroutine), Lemma Sel-DL, Proposition XK-U (exact integer
dependencies at growing degree d ≍ D^{1/5}: Lemmas XK0-U, XK4-U), Lemma C25 (the DL constant equals the factoring
constant), and the cost formula c_DL = max(3a, a+s+max(a,σ), a+s+2^{−2/5}c_U(τ′)).

Base texts: `notes/dlog_results.md` (fixed-degree Theorem DL, Prop XK), `paper2/main.tex` §16, `notes/lv5unif_results.md`,
`notes/recurse_results.md` (Lemma S″, Lemma W, Lemma F), `notes/thmA_proof.md`. Code: `src/dlog25/params.py`.
You may treat Theorem U and its chain as given. Your job is to break what is NEW.

## What to do

Try hard to find a FATAL or MAJOR gap. In particular:
1. **Transfer to N := p.** Does any ingredient of Theorem U use N composite, ω(N) ≥ 2, factor exits, or the square-class law
   of R in a way that breaks for prime N and T uniform on F_p^*? Check select-U → Sel-DL (rounds succeed with probability
   φ(n)/n given the past; size-biased kept field; Wald with random-cost subroutine calls).
2. **Lemma S_p″ and K-U.** Exactness of #{T: E(T)} = #coprime pairs; the window choice; 2Y′² < p.
3. **Hidden labels.** Is every pre-extraction computation really a function of 𝒯_j? Can the randomized shape test
   (FACTOR_U, whose running time and internal recursion depend on its coins) leak label information, e.g. through
   time budgets, alternation with exhaustive search, or restarts?
4. **XK at growing degree.** Re-derive: conjugate bounds [p^{−d}, p]; Blanksby–Montgomery gap and the factor d in
   log M(z) ≤ d‖Log z‖_∞; the coordinate bound (3/2)^{r}2^{r/2}‖v‖/β_d for LLL-reduced bases; Cramer bound R; κ; the
   precision ledger; root approximation; local maximal orders at O(log p) bad primes; Hadamard bound for the rank-profile
   prime; d right-hand sides; torsion with d odd. Find any hidden exponential-in-d or (log p)^{d}-type cost. Is the total
   L^{3a+o(1)} in the L_p[2/5] sense?
5. **Cost and Lemma C25.** Re-derive (6.1) and Lemma C25 (i)–(iii) exactly; check the numbers table with
   `src/dlog25/params.py`; is the claim "3a never binds" justified?
6. **Scope claims.** §7 (no cheap per-log split), §8 (medium characteristic: "does not transfer as it stands"), literature
   (§1): overclaims?

## Output

Markdown. One-line verdict: ACCEPT / ACCEPT WITH FIXES / MAJOR GAP / FALSE. Then an issue table (FATAL / MAJOR / MINOR,
location, why, fix), the weakest link, and answers to items 1–6. Quote lines; give counterexamples or exact estimates.
