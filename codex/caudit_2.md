# Adversarial review of a factoring-barrier COMPLETENESS matrix

A meta-audit claims the following about a factoring-barrier catalogue. Your job: **attack it.** Find (a) any
closure the audit MISSED — a route it lists as "uncovered/open" that is in fact closed by an existing argument or
is trivially factoring/DL-equivalent; and (b) any FALSE GAP — a route the audit calls "illusory/endpoint-
dominated" that is actually a live, genuinely-uncovered path to sub-L[1/3] factoring. Be concrete and skeptical.

Notation: L_N[α,c]=exp((c+o(1))(log N)^α(log log N)^{1−α}); balanced semiprime N=pq, p<q≤2p.

## The audit's claims

VERDICT: the catalogue is a patchwork of ~40 model theorems + ~1000 flagged "open/not covered" notes; it is not
literally complete, but every uncovered corner is (1) endpoint-dominated (only improves the L[1/3] constant, e.g.
Prop SI→1.2114, MNFS/factory→1.6386; or costs N^{1/4} to cash in because every leak gives O(1) factor bits and
Coppersmith needs n/4 bits), (2) secretly factoring/DL-equivalent, or (3) a real coverage gap with no candidate
algorithm. No credible route beats L[1/3]. P(sub-L[1/3]) < 1%.

KEY STRUCTURAL HOLE claimed: the union of single-family barriers is an UNPROVEN EXTRAPOLATION. Leak barriers
(decision trees over Z/N-algebras: PT, NE1, PT-J, AC0-ε) and relation barriers (size-entropy law SS in "model S";
heuristic NFS; listing LS) are proved in DISJOINT models. No theorem bounds a HYBRID that uses a leak to seed a
relation method. leakmap itself says unions are uncovered and integer representatives make any union complete.

EXACT CHECKS the audit relies on:
- Size-entropy law: α=a/(2−b) for candidate size-law λ=A x^a(log x)^{1−a−b}C^b. Sub-1/3 ⟺ 3a<2−b. Realized
  methods sit at α≥1/3 (NFS at a=1/2,b=1/2→1/3). BUT a stylised rank-r law (a=1/2,b=0) formally gives 1/4 and is
  unrealized only heuristically (resbarrier G_r). So SS is a REDUCTION, not a prohibition.
- MNFS/factory are exponent-1/3, constant-only: (64/9)^{1/3}=1.9230, ((92+26√13)/27)^{1/3}=1.9019, 1.6386.
- Prop SI = (16/9)^{1/3}=1.2114, exponent 1/3.
- Coppersmith: ⌈n/4⌉ known top bits of p factor N; O(1)-bit leaks are exponential to exploit.

THE 5 MOST-CREDIBLE UNCOVERED ROUTES the audit lists (most credible first):
1. Cross-family hybrid: a leak biases candidate generation so the effective size-law has a<1/2 (escapes H_S's
   no-per-component-filtering clause). Audit says illusory because known leaks are O(1)-bit, but admits it is not
   closed by a theorem. [This is the one the audit most fears.]
2. Adaptive/non-oblivious smooth-number constructor near kN (beats the oblivious L[2/3]/output floor). Audit:
   constant-only (1.2114) AND DL-equivalent (navigate Theorem DLB).
3. Many-factor-base / super-constant-k sampling outside model S. Audit: covered by heurnfs, exponent 1/3.
4. PRIM nested/retained N-exponent gates (a^N+1)^N. Audit: O(1)-bit leak, cashing in is N^{1/4}.
5. Single modular coefficient a_f(N) of a fixed non-CM cusp form. Audit: open "new view", not a route; native
   cost conjectured ≥ N^{1/2}.

## Your task

1. Is any "uncovered/open" route above actually CLOSED by a standard argument, or trivially equivalent to
   factoring/DL? Name it and the argument.
2. Is any "illusory/endpoint-dominated" demotion WRONG — is there a concrete mechanism (even sketchy) in routes
   1–5, or a route the audit omitted entirely, that could reach exponent < 1/3? Push hardest on #1 (cross-family)
   and on anything the size-entropy law genuinely does not constrain (non-sampling methods; algorithms that use N
   only through an auxiliary object of size > N; adaptive constructions).
3. Does the union of single-family barriers in fact imply a cross-family barrier, by some argument the audit
   missed? Or confirm it is a real hole.
4. Name any route family the audit's 12-row matrix omitted entirely.

Web search allowed. Give your strongest 3 points, most serious first. Do not pad.
