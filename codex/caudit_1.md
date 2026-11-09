# Blind completeness audit of a factoring barrier catalogue

You are an adversarial referee. A research project claims that NO classical algorithm can factor a general
balanced semiprime N = pq below L_N[1/3] (heuristic) / below L[2/5] (rigorous) / below N^{1/5} (deterministic),
and supports this with a large catalogue of "barrier" theorems. Each barrier is a MODEL theorem: it rules out a
class of algorithms within a restricted, formally specified model. Your job: **name every route to sub-L[1/3]
factoring that this catalogue does NOT actually close, and explain precisely why each named route escapes the
stated models.** Be concrete and skeptical. We are hunting for a genuinely uncovered paradigm.

Notation: L_N[α,c] = exp((c+o(1))(log N)^α (log log N)^{1−α}). Balanced semiprime p<q≤2p.

## The catalogue (one line each; these are the model theorems claimed proved in the restricted model named)

RELATION / SIEVE FAMILY
- Size–entropy law (Theorem SS): a single-family SAMPLING method that generates candidates sharing one factor base,
  each candidate tested in full on O(1) smoothness integers, costs exp(C*) where C*²≥2λ(C*)logλ(C*), λ=candidate
  log-size law. Gives α=a/(2−b) for law λ=A x^a (log x)^{1−a−b} C^b. Reproduces Dixon, QS, NFS, SNFS exactly.
  Sub-1/3 needs candidate log-size growing slower in N (a<1/2) or in entropy (b<1/2). MNFS lies OUTSIDE model S
  (shared rational side, V fields on survivors), which is why it reaches 1.902.
- Listing (Theorem LS): if smooth objects are LISTED not sampled, cost is per-output. Doubly-smooth lister ⇒ poly
  time but is factoring-equivalent / gives prime-field DL (Prop DLT). One-sided Bohr-set lister similarly.
- Prop SI: constructing an L[1/3]-smooth integer in EVERY short interval [kN,kN+L[2/3,ν]] ⇒ L[1/3,1.2114]. Open;
  a one-sided Schnorr-CVP reformulation; no sub-sieving constructor known (every oblivious constructor costs L[2/3]
  per output, Theorem J).
- Heuristic NFS: nothing beats MNFS 1.902 / exponent 1/3 across polynomial selection, lattice sieving dimension,
  output-sensitive smooth-pair enumeration (validated on CADO-NFS).

COVERAGE / DESIGNS (deterministic, exponential regime)
- N^{1/6} true-factor floor for local collision designs (T1/T2, multimod Theorem A); N^{1/5} with one modulus
  (cycloharvey LB); coverage cost ≥ N^{1/5}(λ⁴/J)^{1/5}. Worst case < N^{1/5} needs coverage exponent c>3/8
  (dchworst Theorem W), shown to need a 7/25 saving that phase-discarding methods cap at 4/25 (cov38b Theorem B).

DESCENT LADDERS (rigorous randomized)
- One-level Dixon descent of a random square ⇒ L[2/5] (Theorem U). Two-level ⇒ L[3/8] conditional/then via
  published inputs. j levels ⇒ (j+1)/(3j+2) → 1/3 as j→∞, but depth ≥4 fails: required precision falls below the
  natural covariance. Rigorous L[1/3] needs ERH + far stronger all-target relation supply (nested Theorem NR).

FLOOR-SUM / ANALYTIC / CIRCUIT
- Floor-sum ladder at degree j≥3 = Coppersmith known-bits barrier (topblock). Analytic sums, short circuits
  (slpcirc Theorem SP), period-finding (Theorem T: only large-order multiplicative-character observables peak,
  which leads back to smoothness). Pochhammer/factorial SLPs ⇒ factoring in P/poly but contradict τ-conjecture.

LEAKS WITHOUT SMOOTHNESS (where does "information cost splitting" stop?)
- Every known "leak" (way to make p,q behave differently without splitting N) is one of three families:
  group-order/smoothness (M), archimedean (A), or quadratic genus (Q). Proved for: N-dependent exponent trees
  (NE1), AC⁰ on bits of N (target ε), pure trees + Jacobi symbols (PT-J). Integer representatives are
  computationally complete (no PT-type theorem). Strongest OPEN model "PRIM": uniform programs over fixed Z/N
  algebras whose only N-dependence is retained/nested exponent gates x↦x^{E(N)} from a fixed menu, + gcd/zero
  tests + Jacobi of computed values. One-layer (NE1) and terminal-Jacobi parts proved; nested/retained gates OPEN.

REPRESENTATION / SELF-REDUCTION / NATIVE
- Addition in prime-exponent representation ⇔ factoring. No modulus-changing self-reduction in algebraic classes
  (SR1,SR2). 55 factoring-equivalent invariants (φ,σ,class numbers,modular coeffs,K-groups,point counts): none has
  a native algorithm beating N^{1/5}; open view: one coefficient a_f(N) of a fixed non-CM cusp form.
- #E(Z/N) for a fixed curve factors N. Deterministic subexp under ERH reduces to one constructive lemma DSR; every
  black-box derandomiser fails on it.

CONNECTIONS: ~300 surveyed across 40 fields, none a route; seven failure mechanisms (restatement, aggregate mod m,
ring-arithmetic-mod-N, group-order, sampling-precision 1/p, CRT tensor, distributional-only).

## Your task

1. List routes/paradigms to sub-L[1/3] (or sub-2/5, sub-3/8) factoring that are NOT closed by ANY model above.
   For each: (a) the mechanism; (b) which model it would have to live in and why it escapes that model's
   hypotheses; (c) a rough exponent estimate; (d) whether it is secretly factoring-equivalent or secretly
   endpoint-dominated (i.e. the gap is illusory).
2. Pay special attention to: NON-OBLIVIOUS / ADAPTIVE constructions; mechanisms that use N only through an
   auxiliary object of size > N; HYBRID routes combining two individually-barred families; anything the
   size–entropy law does NOT constrain (it only covers single-family sampling with O(1) smoothness integers and
   one shared factor base — what about many factor bases, super-constant k, non-sampling?); whether the union of
   single-family barriers really implies a cross-family barrier.
3. Give your 3 best candidate uncovered routes, most-credible first, with the cheapest check that would settle
   plausibility.

Use web search if helpful. Be concrete; prefer a precise escape over a vague one. Do not pad.
