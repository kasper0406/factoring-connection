# floor8: fast mod-8 evaluation of a 2D floor sum with a quadratic congruence N | Q_0(x,y)? (2026-11-26)

**Question (user).** For a fixed binary quadratic form Q_0 (disc D) and N = pq, can the 2D count
S_B(N) = #{(x,y) in B : N | Q_0(x,y)}, or a genuine floor variant S'_B = sum floor(x/2) over the same points, be
computed mod 8 in poly(log N) WITHOUT splitting the congruence (no factoring, no sqrt(D) mod N)?

**Process.** GPT-6.2-astra via codex (session 01a218e5-a20b-7ef2-b395-c7df38b022c3). Brief: `codex/floor8.md`. Output:
`codex/floor8_codex.out`. Artifacts (incl. a verified factoring routine): `experiments/floor8/`.

## Verdict: no useful fast mod-8 evaluation — a variable-size or floor mod-8 evaluator FACTORS N (proved + verified)

### Fast (proved / known)
- **Complete-period count is fast mod 1, 2, 4 for every N** by binary quadratic Gauss-sum reciprocity
  (Deloup–Turaev / Landsberg–Schaar; |D| terms, needs only N mod 8|D|, Jacobi (D|N), gcd — no factoring). The density
  r(N) = #{Q_0 ≡ 0 mod N} is multiplicative:
  r(p^e) = p^{e-1}((e+1)p − e) if (D/p)=1, else p^{2⌊e/2⌋}; both ≡ 1 mod 4.
- **mod 8 is fast for D = −4** (x²+y²): every good split prime is 1 mod 4.
- The Weil-index phase and genus character are the fast computable shadows (H1 confirmed).

### Hard (proved; contains the u_2 bit)
- For D = −23, complete-period mod 8 is
  r(N) ≡ 1 + 4·#{p | N : p ≡ 3 mod 4, (D/p) = 1} (mod 8),
  which on the (D/N)=1, N≡3 mod 4 stratum is exactly u_2(N) = [D is a square mod N] = [R(N) = 4] (kern/coef). So even
  complete-period mod 8 holds the factoring-hard bit. (H2 confirmed, via the two local characters, not the class-group
  2-part directly — for D=−23 the class group is C_3 with trivial 2-part; the bit is principal-vs-nonprincipal split.)

### Variable-size / floor: a fast evaluator IS a factoring algorithm (proved + verified)
- **Reflection identity:** T(k) := 2 Σ_{i≤k} h(i) = r − 1 + S(k+1) − S(N−k), all queries X ≤ N.
- T(k)/4 mod 2 = parity of ⌊k/p⌋+⌊k/q⌋ (both split) / ⌊k/p⌋ (one inert) / 0 (both inert). At k = ⌊N/2^j⌋ these are
  the bits of p XOR q (resp. q). **Bisecting on a differing-parity endpoint lands the final jump on a multiple of a
  factor → gcd recovers it.**
- **Genuine floor ⌊x/2⌋ (H3):** exact second difference E(k+1) − E(1) = 3N·T(k); 3N odd ⇒ W mod 8 recovers T mod 8 ⇒
  same factoring reduction. The floor does NOT escape.
- Even variable-square **parity** (mod 2) factors a subfamily (p inert, q split, (6/q)=−1), via Gauss's lemma on two
  lines mod q.

### Measured (pre-registered, 40–80 bits, D = −23 and −4; ground truth from factors)
- Both complete-period identities passed every case (incl. W(N) = (r−1)(N−1)/4).
- Registered public features (J = (D/N) mod 8, K = 12N·s(2,N) mod 8, N mod 8) gave **50% balanced accuracy** on
  R = 4 vs 8; the planted-signal control rose to 100%. r mod 8 does not determine R(N) and vice versa.
- **Factoring routines verified:** registered bisection 365/365 and 355/355 split/split; folded bit-query 1,116/1,116
  and 1,140/1,140 of semiprimes with any split prime; fixed-floor second difference passed all 3,072 held-out pairs;
  variable-parity subfamily 364/364.

## Dedekind reciprocity (why the floor doesn't rescue it)
Classical s(h,k) reciprocity evaluates a *supplied linear slope* by Euclid. Dedekind–Rademacher / cotangent
generalizations take their linear arguments as input; applied to the congruence lattice they need the splitting first.
The 2-part of s(h,k)'s reduced denominator is small (12k·s ≡ k+1−2(h/k) mod 8 for odd k), so mod-8 control is possible
in principle, but evaluating the unsplit quadratic-congruence sum is exactly the factoring reduction above.

## Sharpest open sub-statement
A SINGLE query at the prescribed short square X ≈ √N for D = −23 (unit weight or fixed ⌊x/2⌋): the adaptive reductions
all use many X and leave the one-shot short-square count mod 8 open. For x²+y², that short-square count is exactly
1 + 4·[both primes split]. Astra's probability of a fast factor-free mod-8 algorithm uniform in X: ~1%.
