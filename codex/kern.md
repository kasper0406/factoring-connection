# kern: Frobenius-free counts on the ℓ-torsion algebra over Z/N, and the OR-structure of zero-detection

Context (this session's coef/dig work): A = (Z/N)[x]/(ψ_ℓ(x)) (or the full ℓ-torsion algebra with the curve
equation) is computable; Frobenius on it is not (producing componentwise Frobenius factors N — Prop FROB,
notes/nonsmooth_results.md:313). The open gap is computing a NUMBER (u_ℓ(N) = det(1−g_p)det(1−g_q) mod ℓ, or
whether ℓ | M) without producing the operator. Repo (do NOT edit; scratch in your cwd): notes/coef_results.md,
notes/dig_results.md, notes/tracedata_results.md (Theorem EXT, Prop NA/IL, mixed points), notes/ecount_results.md,
notes/nonsmooth_results.md, notes/prim_results.md, notes/factorbits_results.md (Theorem PT).
Two focused ideas from the user:
K1. The multiplicative map μ_N: a ↦ a^N is a group endomorphism of A^× (not a ring map). On A_p^× ≅ ∏ F_{p^{k_i}}^×
    (k_i = orbit lengths of Frobenius on the roots of ψ_ℓ mod p — these orbit lengths ENCODE the Frobenius class g_p up
    to the information in its action on lines/points), its kernel has size ∏ gcd(N, p^{k_i} − 1), similarly on the q
    side. Questions: (a) exactly what Frobenius data (orbit structure of g_p, g_q on P(E[ℓ]) or E[ℓ]∖0) is determined by
    the multiset of k_i, and does that determine u_ℓ? (b) Is |ker μ_N| (or |ker μ_m| for chosen m, or the order of
    A^×, or the number of solutions of a^m = 1, or |A^× / (A^×)^m|) computable from N without factoring — e.g. via
    counting arguments, Jacobi-type symbols in A (the "Jacobi symbol" of an element of a finite étale Z/N-algebra
    w.r.t. its norm/discriminant), resultants, or random sampling (estimate |A^×[m]| by sampling a ↦ a^{|A^×|/m}...
    but |A^×| is unknown — is it?). Note |A^×| = ∏ (p^{k_i} − 1)(q^{k_j} − 1) is itself a "group order" leak (family M);
    knowing it would factor N? Decide. (c) Any sampling-based estimator of the k_i statistics (e.g. frequency with
    which a^{N^j} = a for random a in A, for j = 1..ℓ²; or the period of j ↦ a^{N^j}) — does a^{N^j} relate to Frobenius
    powers on each side? (On the p side a^{N^j} = a^{p^j q^j}; when does it equal a? Give exact conditions in terms of
    k_i, p, q.) Run a cheap pre-registered experiment (40–80 bits, ground truth from factors) for the most promising
    statistic, with gcd-split flagging.
K2. Zero-detection u_ℓ(N) ≡ 0 ⟺ ℓ | M ⟺ E(F_p)[ℓ] ≠ 0 OR E(F_q)[ℓ] ≠ 0. Z/N-valued torsion points see the AND (both
    components); the OR needs mixed CRT points, which are factoring witnesses (EXT). At ℓ = 2 for y² = x³−x+1 this is
    exactly "R(N) = #{s²+st+6t² = N} ∈ {0, 8} vs 4". Questions: (a) Is there ANY computable test that detects the OR
    without constructing a mixed point (e.g. counting arguments over A, the number of idempotents, |A^×| parity,
    discriminant/Jacobi of the norm form, trace forms of A: the trace form Tr_{A/(Z/N)}(xy) has a discriminant and a
    Hasse–Witt / Hasse invariant computable mod N — does its square class or a related quadratic invariant detect
    splitting types? (Serre's trace-form formula relates the trace form's invariants to Stiefel–Whitney classes of the
    Galois representation — a non-abelian invariant with abelian-looking values!) (b) For ℓ = 2 specifically, does the
    trace form of A = (Z/N)[x]/(x³−x+1), or of the 2-torsion algebra, or the Hasse invariant of some explicit quadratic
    form built from N and the cubic, compute (or correlate with) R(N) = 4? Run an exhaustive / pre-registered check.
K3. Verdict: what is proved, what is refuted, the sharpest remaining open statement, probability.
≤ 2000 words; label proved / known / heuristic / measured / speculative.
