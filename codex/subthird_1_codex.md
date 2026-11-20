(Common preamble: read /Volumes/git/factoring-connection/codex/subthird_common.md first and follow it.)

## Your focus (consultant C): constructive attack on the G_r / type-D crack, with experiments

1. Restate the grdeg open case exactly: fixed small subfield k (degree e), K/k of degree n, linked O_k-module of
   relative rank m ≥ 3 with non-power generators, and the L[1/4] ledger (grdeg §0 deliverable #1; midk §0). Identify
   exactly which counting statement (G_r independence / H5_r) has to fail and in which direction for α < 1/3.
2. Try hard to CONSTRUCT such a family: polynomial selection for N giving m ≥ 3 linked generators over k with cheap
   linking (no factoring of N needed), small heights, and many independent low-norm box points not confined to
   lower-rank sublattices or killed by unit action. Use the repo's code (src/grdeg, src/rank3, src/rankdegen) as a
   base; run toy experiments (N up to ~60–100 bits; exact arithmetic; Python/PARI/Sage if available) measuring the
   actual exponent of independent relations vs the G_r prediction.
3. Either exhibit a family whose measured relation yield beats the G_r/L[1/3] prediction with a scaling argument, or
   prove/argue the obstruction (e.g. a relative R1 over O_k with error o(x^{3/4}) — note a Bost-slope version with
   O(dn) error is vacuous at the L[1/4] balance).
4. Separately, give ≤ 2 other mechanisms you think have a real (>0.5%) chance at α < 1/3, with ledgers.
