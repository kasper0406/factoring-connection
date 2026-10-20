# cadobench round 1: practical NFS speedup levers 2020–2026, and what CADO-NFS already does

You are a literature and engineering expert on the general number field sieve (GNFS) as it is actually
implemented. Use web search. Cite sources (paper, ePrint number, CADO-NFS commit or file, mailing-list post)
for every claim. Where you are unsure, say UNVERIFIED.

Context. We are building a real-size testbed around CADO-NFS (current git master, Oct 2026) on an Apple M4
laptop: 4 performance + 6 efficiency cores, 16 GB unified memory, NEON, SME in streaming mode only, Metal GPU.
We factor balanced semiprimes of 60–110 digits and profile `las`. We want *practical, measurable,
constant-factor* speedups. We do not want to redo known work. Our heuristic analysis already says that, at the
L[1/3] level, cheaper smoothness tests cannot move the MNFS/GNFS constant, and only output-sensitive
enumeration of smooth pairs or sub-quadratic sparse linear algebra would. So this round is about engineering
constants.

Questions. Answer each in a short section.

1. **The practical levers, 2020–2026.** List the levers that measurably sped up real GNFS records or
   implementations in this period, with rough reported gains. Include at least:
   - polynomial selection (Kleinjung 2nd/3rd generation, size-optimisation and root-optimisation, Gaudry–
     Thomé–Zimmermann et al., the RSA-240/250 record polys);
   - lattice sieving (bucket sieve, small-sieve, "las" vs Franke–Kleinjung lasieve4 vs others; 3 large primes on
     one side (3LP); sieving both sides vs one; I/A parameters);
   - cofactorization (ECM/P−1/P+1 chains, Kruppa's strategies, "cofactorization strategies" optimisation by
     Bouvier–Imbert, batch smoothness / Bernstein product trees, CADO `batch=1`);
   - filtering (singleton/clique removal, merge improvements by Bouvier, Cavallar, Zimmermann et al. 2021+);
   - linear algebra (block Wiedemann parameters, `bwc` with MPI, GPU block Lanczos in msieve-gpu, block size,
     the RSA-250 / DLP-240 matrices);
   - square root.
2. **What CADO-NFS already implements** for each lever above, and which parameters control it (e.g. `lpb0`,
   `lpb1`, `mfb0`, `mfb1`, `lambda`, `ncurves0/1`, `batch`, `batchlpb`, `batchmfb`, `I`/`A`, `qmin`, `rels_wanted`,
   `tasks.sieve.adjust_strategy`, `bwc` `interval`, `m`, `n`). Which parts are known to be weak or untuned?
3. **GPU and SIMD offload.** What is known (2020–2026) about GPU NFS sieving or cofactorization: GPU-ECM
   (Bos–Kleinjung, Miele et al., "cuda-ecm", ECM on GPUs 2020–2025), GPU lattice sieving, GPU block Lanczos/
   Wiedemann? Measured gains on real factorizations? Anything Apple Silicon/Metal-specific? Any ARM NEON/SVE
   work inside CADO (e.g. bucket sieve, `fill_in_buckets`, norm initialisation)?
4. **Batch smoothness in practice.** When does CADO's `batch=1` path (product trees over large-prime ranges)
   beat per-candidate ECM chains? Reported numbers? What were the RSA-240/250 and DLP-240 settings?
5. **Parameter tuning at c90–c110.** Are CADO's default parameter files (`parameters/factor/params.c90` …
   `params.c110`) known to be near-optimal, or have people reported 10–30% gains by retuning lpb/mfb/lambda/I?
   Cite forum posts (mersenneforum) if that is where it is.
6. **Ranked recommendation.** For a single M4 laptop at c90–c110, rank the 5 levers most likely to give a
   measurable ≥5% wall-clock gain over CADO defaults that we could test in a week, with your estimate of the gain
   and the evidence for it. Separately, flag levers that are already exhausted (CADO does them well).

Keep the answer under 2500 words. Be concrete, give numbers, and separate KNOWN (cited) from your own estimates.
