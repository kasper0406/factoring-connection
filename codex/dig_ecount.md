# dig A: a native algorithm for #E(Z/N) of ONE fixed curve (Schoof at composite N)

You ranked this the most likely home of a poly-time factorer (notes/inp_results.md). Dig in. Repo (do NOT edit;
scratch in your cwd): notes/tracedata_results.md (Theorem EXT: one curve's count M = (p+1−a_p)(q+1−a_q) factors every
good N), notes/ecount_results.md (Okamoto–Uchiyama; one count ⇒ factoring), notes/nativemod_results.md
(Edixhoven–Couveignes at composite N is FROB/PRIM; Remark 2.2.4 prime/composite dichotomy), notes/nonsmooth_results.md
(Prop FROB: x ↦ x^N on Z/N-algebras is not Frobenius), notes/prim_results.md, notes/factorbits_results.md (Thm PT),
notes/native_results.md, notes/geometry_results.md (X_0(N) point counts), notes/conn_alg.md.
Questions (≤ 2000 words; label proved / known / heuristic / speculative):
E1. Exactly where Schoof (and SEA, Pila's generalization, Edixhoven–Couveignes, p-adic/Kedlaya/Satoh-type and
    deformation methods, Harvey's average polynomial time, Lauder's deformation, AGM) uses that the base ring is a field
    or that p is known: the Frobenius endomorphism π on E[ℓ], its characteristic polynomial, x^p mod (ψ_ℓ, p). For each
    method, what replaces the missing ingredient over Z/N, and why the obvious substitute (x^N, CRT-componentwise
    Frobenius, Witt vectors, deformation from a known curve) fails or is factoring-equivalent.
E2. Reformulations that do not need Frobenius: M mod ℓ via the action of SOMETHING computable mod N on E[ℓ] (e.g.
    the Weil pairing, isogeny volcano structure, modular polynomials Φ_ℓ(j, X) mod N, division-polynomial resultants,
    the ℓ-adic Tate module via lifting to Z/N^k, L-function / modular-form coefficient approach a_N vs the count
    (note: M is NOT a coefficient — M = N + 1 − a_p q − a_q p + a_p a_q ... write it in terms of computable quantities)),
    CM curves (where #E(F_p) is a norm form in the Frobenius π ∈ O_K, so M = N(π_p − 1)N(π_q − 1) — computing it ≡
    representing primes by forms), supersingular/anomalous special curves, and higher-dimensional analogues (Jacobians,
    abelian surfaces, Kummer). For each: the exact obstruction or a precise open subproblem.
E3. Partial information: is M mod 2, M mod 3, M mod small ℓ computable from N without factoring (e.g. from the
    structure of E[2](Z/N) = roots of the cubic mod N, E[3] via 3-division polynomial roots mod N)? Note M mod ℓ for
    small ℓ alone gives little; but how much does Theorem EXT need — the full M, or M mod a modulus > N^{1/4}·something?
    Precisely quantify the recovery side's requirements.
E4. Run a cheap pre-registered experiment: e.g. test whether M mod small ℓ is predictable from N-computable features
    (root counts of division polynomials mod N, Jacobi symbols, etc.) beyond what CRT-structure trivially forces, at
    40–80 bits, with exact counts from the known factors as ground truth. Report.
E5. Verdict: the single most precise open subproblem whose solution would give native #E(Z/N), and your probability.
