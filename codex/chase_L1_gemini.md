(Read /Volumes/git/factoring-connection/codex/chase_common.md.) LEMMA 1 (H / branch selection).
Context: experiments/hatch/hatch_h1c/report.md (codex's quadratic rational reconstruction: once a divisor/branch is
chosen, witnesses are recovered fast; 5.96 vs 10.27 ms per witness at 96 bits; the open step is selecting the branch
BEFORE enumerating L[1−β] divisor products), notes/subthird_results.md (M2 class-group 2-D CRT constructor),
notes/pilist_results.md, notes/crtdecode_results.md, notes/sihard_results.md (Thm J), notes/hatch_results.md.
Target: an algorithm that, in poly(π(B)) (or L[β]) work, selects a divisor branch (a product of split-prime roots /
a CRT branch) whose quadratic reconstruction yields a B-smooth witness, beating the Dickman rate — or a proof that
branch selection is as hard as the original search (e.g. reduce random-branch selection to Theorem J or to a
counting bound). Consider: meet-in-the-middle over root branches, lattice enumeration with pruning keyed to partial
norms, sieving-in-the-branch-tree (does partial-product information predict full smoothness?), sampling branches
proportional to a computable proxy. Measure on toy sizes (python3) whether any proxy beats random branch choice at
fixed work. Use only read-only shell commands plus python3 in your cwd. Final answer ≤ 2000 words, plain text.
