(Read /Volumes/git/factoring-connection/codex/guide_common.md.) Role: BUILD AND RUN.
1. Pipeline: for several N (60–110 bits, balanced semiprimes, ≥ 6 moduli, split train/test by N), get NFS polynomials
   (CADO polyselect or base-m), enumerate a large set of candidates (≥ 10^6 labeled points total; both full-box (a,b)
   and special-q lattice regions), compute both norms, strip primes ≤ y, label B-smoothness of the cofactors exactly
   (factor with PARI or ECM up to the needed size; or use las output plus independent verification).
2. Baseline: calibrated logistic/isotonic on log-cofactor sizes (+ the y-local information). Test models: gradient
   boosting and a small MLP on leakage-safe features (see preamble), plus a "history" feature set (outcomes of other
   candidates, excluding same-residue-class neighbours for primes in (y, B]).
3. Pre-register (file in cwd, before running): metrics, held-out protocol, thresholds (e.g. "signal" = held-out
   log-loss improvement > 3σ over baseline and replicated on a second N-family). Run. Report exactly, with leakage
   checks (permutation tests; shuffled-label controls; a planted-signal positive control to show the pipeline CAN
   detect a real correlation of the size that would matter).
4. Ledger: if any signal is found, estimate how much yield gain per unit work it gives and whether it could compound
   to an exponent change (it needs to grow with u, not stay O(1)).
Keep it within ~2 hours of compute. Write everything to files in your cwd; final summary ≤ 1500 words.
