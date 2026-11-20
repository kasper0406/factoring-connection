# lowerbound round 2: cross-review and merge (rules: codex/lowerbound_common.md)

Round-1 (paths relative to /Volumes/git/factoring-connection):
- C (codex): codex/lowerbound_1_codex.out and full proof codex/lowerbound_1_codex_proof.md (supermartingale
  composition theorem; size-profile hypothesis (2) λ(k) ≥ A√(x log(k+1)) ⇒ L[1/3,(4A/3)^{2/3}]; proved profile inputs
  for unions of binary forms (Lemma C2) and a KH subclass; GNFS/MNFS tightness).
- G (Gemini): codex/lowerbound_1_gemini.out (literature) and codex/lowerbound_1_gemini2.out (model M_G with Axiom NP,
  SM1/SM2/H_sel, LP_exc^rel, Axiom R; Theorem G; Step 2 removing Assumption T via analyst cutoff B_* = T²(log N)^C;
  batch cap K ≤ N^{1/2−ε} justified by adversary A1; adversaries A2–A6).

You are one of the two. ≤ 1500 words:
Z1. Check the other's central proof step(s) hard. For C: the supermartingale/Markov bound (1), the dyadic size-entropy
    consequence and constant (4A/3)^{2/3}, the sharing/credit rules (can they be gamed?), and whether hypothesis (2) is
    really implied for GNFS/MNFS by the proved inputs or merely matches the heuristic profile. For G: Step 2 (does the
    SM2 birthday bound really force every large prime of a dependency into a forced lattice modulus, and does the last-
    step argument then give a fully B_*-smooth unforced cofactor? what about large primes created by forced moduli of
    earlier steps?), the universal constant (1/18)^{1/3}, the special-N floors 1.2114/1.5263, the batch-cap argument
    (union bound T·K ≤ N^{1−ε}), the Axiom NP formulation, and the ω = 0 constants 1.8595/1.6363.
Z2. Is adversary A1 (dense-batch factory with K ≈ N/log N giving polylog amortized cost) correct? Is K ≤ N^{1/2−ε} the
    right cap, or is a smaller/larger cap needed?
Z3. Merge: write ONE theorem statement (model + assumptions + conclusion) that takes the best of both — C's adaptive
    supermartingale composition with G's N-publicity axiom, assumption-T removal (if it survives) and batch cap — listing
    exactly which assumptions are proved for which algorithm classes and which are heuristic/model assumptions
    (SM/H, R, profile (2)/LP_exc^rel, H_sel). Keep it tight enough that a referee could check it.
Z4. The 3 points a hostile referee would attack first, and whether each is fixable.
