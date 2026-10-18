# genus2rh: RH-conditional rigorous L_p[1/2, 2] factoring via genus-2 Jacobians

Agent `genus2rh`, phase 6, 2026-10-18. Source idea: codex/r18_ideas.out §3.
Files: src/genus2rh/{cantor,hecm,diagnostic,experiment_g1}.py, tests/test_genus2rh.py,
notes/genus2rh_prereg.md, experiments/raw/genus2rh_g1{,b}.jsonl.

## 0. Verdict

**V1 (proved, conditional on RH for ζ; uses one refereed theorem and one unrefereed preprint).**
The statement is correct. Theorem R (§4): assume RH. A Las Vegas algorithm, given a composite
N that is not a prime power, outputs a proper divisor of N in expected time

  L_p[1/2, 2 + o(1)] · (log N)^{O(1)},  where p is the least prime factor of N.

The o(1) is O((log log p)^{-1/2}) as p → ∞, uniformly in N. With schoolbook arithmetic the
polylog factor is (log N)^{5+o(1)}. The proof combines:
- Lenstra–Pila–Pomerance II, Thm 1.1 (refereed, Proc. LMS 2002);
- Younis, arXiv:2409.05761v1, Thm 1.3 (RH; unrefereed, v1 only, as of 2026-10-18);
- Canfield–Erdős–Pomerance 1983;
- a genus-2 transcription of Lenstra's 1987 ECM separation argument (Props A, B below).

Every loss in the accounting is L_p[1/2, o(1)] or (log N)^{O(1)}. No step costs a power of p
(§5, audit table).

**V2 (proved; elementary corollary of V1).** Suppose p ≤ N^β. Then the expected time is
L_N[1/2, 2√β + o(1)], uniformly as N → ∞. Running the Lenstra–Pomerance 1992 algorithm in
parallel gives L_N[1/2, min(1, 2√β) + o(1)] for every N. The constant is below 1 exactly when
β < 1/4. Example: 2/√5 ≈ 0.894 at β = 1/5.

**V3 (novelty: the mechanism and the constant are NOT new; the RH statement as a theorem is
not found by us).**
- (a) The constant 2 is LPP's own. LPP I (Phil. Trans. 1993, p. 398) states that "under a
  similar hypothesis [on smooth numbers in short intervals], the expected run time of the
  hyperelliptic curve method, with optimal choice of parameters, is actually at most
  L_p[1/2, 2+o(1)](log n)^2".
- (b) Younis 2024 is, under RH, a hypothesis of exactly the kind LPP needed. We need only its
  lower-bound half, at θ = 3/4.
- (c) LPP's algorithmic analysis ("part III") was never written. Pomerance (IMA Abel
  Conference talk, Jan 2011, slide 8) calls it '"to appear" (actually, "to be written")'. He
  says the method "stands (or will stand when completed)". §3 below supplies that accounting.
  It is a routine adaptation of Lenstra 1987 §2, with one genus-2 twist: LPP II Thm 1.1 is
  applied at the second prime q for the separation condition.
- (d) We found no paper that states Theorem R, or any RH-conditional rigorous L_p[1/2, c]
  factoring bound. Our searches covered citations of Younis, LPP, Pomerance 1987, Lenstra 1987,
  Lenstra–Pomerance 1992, Hafner–McCurley and Granville's survey (§7). By rule 5 we claim only
  "no prior formulation found by us".
- Assessment: this is a straightforward corollary. It might make a short note or a remark,
  but it is not a research result. The L_N corollary is the only part that "beats L_N[1/2,1]",
  and it does so only for inputs with a factor below N^{1/4}, under RH.

**V4 (side result, proved modulo the OCR caveat in §6.3).** Prop A + Prop B with LPP I Thm 1.3
(unconditional) in place of Younis gives an unconditional expected time L_p[2/3, O(1)](log N)^{O(1)}.
This is LPP I's Thm 1.1, whose proof depended on the unwritten part III. §3 thus completes it.

**V5 (optimisation; proved for this analysis model).**
- s = 1 is optimal (cost s + 1/s).
- Genus g ≥ 3 is worse: the cost is √(2g) even if LPP-type fibre bounds existed, and they do
  not exist.
- Genus 1 would give √2. It needs θ = 1/2, which Younis-type RH methods cannot reach.
- Split Jacobians, stage-2 continuations and large-prime variants change only L^{o(1)} factors.
- So 2 is the floor of this route (§6).

**V6 (unconditional / almost-all; not available).** No unconditional theorem known to us
gives L[1/2]-smooth density in the LPP interval, even for almost all p. Matomäki–Radziwiłł
needs fixed u. Younis Thm 1.1 and Jain need y ≥ exp((log x)^{2/3+ε}). See §7.

## 1. Notation

- L_x[a, c] = exp(c (log x)^a (log log x)^{1-a}).
- P^+(n) is the largest prime factor of n, and Ψ(x, y) = #{n ≤ x : P^+(n) ≤ y}.
- I_p = [p² − ½p^{3/2}, p² + ½p^{3/2}].
- λ(t) = (log t)²(log log t)².
- A tuple is τ = (a, b, f_1, …, f_5) ∈ (Z/N)^7. It defines f_0 := b² − Σ_{i≥1} f_i a^i,
  f = Σ f_i X^i, and the point (a, b) on y² = f(x).

## 2. The two analytic inputs: exact statements

**LPP II = H. W. Lenstra Jr., J. Pila, C. Pomerance, "A hyperelliptic smoothness test, II",
Proc. London Math. Soc. (3) 84 (2002) 105–146, Theorem 1.1** (verbatim modulo typesetting; from
math.dartmouth.edu/~carlp/PDF/paper129.pdf, checked against people.maths.ox.ac.uk/pila/HyperII.pdf):

> Let p be a prime number with p > 8100, and let F_p denote the prime field of cardinality p.
> Then for all but at most p integers z in the interval [p² − ½p^{3/2}, p² + ½p^{3/2}] there
> are at least p^{9/2} / (24,000 · (log p)² · (log log p)²) quintic polynomials f ∈ F_p[X] with
> non-vanishing discriminants such that #J_f(F_p) = z.

- Quantifiers: fixed p > 8100; the exceptional set has size ≤ p; the bound holds for every
  non-exceptional z.
- "Quintic" means degree exactly 5, not necessarily monic. There are ~p^6 of them, so the
  average fibre is ~p^{9/2}/8 over the full Weil interval.
- J_f is the Jacobian of the smooth genus-2 curve y² = f(x).
- LPP II states only lower bounds. They write that an upper bound "of the same form" is
  believed but not proved. We need no upper bound.

**Younis = K. Younis, "Asymptotics for smooth numbers in short intervals", arXiv:2409.05761v1
(9 Sep 2024; no later version, no journal reference as of 2026-10-18), Theorem 1.3:**

> Assume the Riemann Hypothesis is true. Let 1/2 < θ ≤ 1. There exists a constant K = K(θ) > 2,
> as well as a function f_θ(x, y) which is independent of h, such that the estimates
>   (Ψ(x+h, y) − Ψ(x, y))/h = α Ψ(x, y)/x · (1 + O_θ(min{1/u + h log(u+1)/(x log y), log(u+1)/log y}))
> and [a second estimate with main term f_θ] hold uniformly for x^θ ≤ h ≤ x and (log x)^K ≤ y ≤ 2x.

- Here u = log x / log y, and α = α(x, y) is the saddle point.
- By Hildebrand–Tenenbaum 1986 (Lemma 2; Younis (1.5)), α = 1 − (log(u log(u+1)) + O(1))/log y
  for log x < y ≤ x. In our range this gives α → 1.
- §7.1 of Younis (the RH proof) derives the version without α, (7.2):
  (Ψ(x+h,y) − Ψ(x,y))/h = Ψ(x,y)/x · (1 + O(log(u+1)/log y)). It holds for
  exp((log log x)^{5/3+1/100}) ≤ y ≤ x^{1/2−ε}. We use either form.
- Quantifiers: θ is fixed, and the implied constants depend on θ only. We take θ = 3/4.

**CEP** (Canfield–Erdős–Pomerance, J. Number Theory 17 (1983), Cor. to Thm 3.1): uniformly for
u → ∞ with y ≥ (log x)^{1+ε}, Ψ(x, y) = x · u^{−u(1+o(1))}. The o(1) is O(log log u / log u).

**Weil.** Let C/F_r be a genus-2 curve with a squarefree quintic model y² = F(X), where r is odd.
Then #C(F_r) = r + 1 − a with |a| ≤ 4√r, and #J(F_r) ∈ [(√r − 1)^4, (√r + 1)^4]. The model
has one point at infinity, so the number of affine points is A(F) = #C(F_r) − 1 ∈ [r − 4√r, r + 4√r].

## 3. The algorithm and its algorithmic accounting

### 3.1 Algorithm C (one curve; operational, src/genus2rh/hecm.py:one_curve)

Input: odd N with no prime factor ≤ 8100; parameters v (a guess for p) and w (smoothness bound).

1. Draw τ uniformly from (Z/N)^7. Set f_0 := b² − Σ_{i=1}^5 f_i a^i, so that f(a) = b².
2. Let g = gcd(f_5, N). If 1 < g < N, output g. If g = N, stop with no output.
3. Pass to the monic model F(X) = f_5^4 f(X/f_5), with point (A, B) = (f_5 a, f_5² b).
   - This is an isomorphism over Z/N because f_5 is a unit.
   - disc F = f_5^{12} disc f.
4. Let g = gcd(disc F mod N, N). If 1 < g < N, output g. If g = N, stop.
5. Set D := (X − A, B), a Mumford pair: u = X − A is monic and u | F − B².
   - Let e(r) := max{m : r^m ≤ (⌊√v⌋ + 2)^4}.
   - For each prime r ≤ w in increasing order, repeat e(r) times: D ← [r]D. This uses
     double-and-add, and each group operation is Cantor's composition-and-reduction over Z/N.
   - Every polynomial is normalised. Coefficients ≡ 0 mod N are stripped, and the top
     coefficient c is tested.
   - If 1 < gcd(c, N) < N, output gcd(c, N). This is the only way a "NonUnit" arises.
6. Otherwise stop with no output.

**Lemma 1 (sampling; proved).**
- (a) For each prime r | N, τ mod r is uniform on F_r^7.
- (b) For distinct primes p, q | N, τ mod p and τ mod q are independent. Reason:
  Z/N → Z/pq ≅ F_p × F_q is a surjective ring map with fibres of equal size.
- (c) Fix a squarefree quintic f over F_r. Exactly A(f) tuples mod r have f_1, …, f_5 equal
  to the coefficients of f and f(a) = b². Their points are exactly the affine points of
  y² = f(x). So Pr[τ mod r gives f] = A(f)/r^7 ≥ (1 − 4r^{−1/2}) r^{−6}.

**Lemma 2 (outputs are proper divisors; proved).**
- Each output is gcd(c, N) with c ≢ 0 mod N and gcd(c, N) ≠ 1, so it is a proper divisor.
- Steps 2 and 4 output only when 1 < g < N.

**Lemma 3 (reduction lemma for Cantor over Z/N; proved).** Suppose a run of steps 3–5, or any
prefix of step 5, ends without a NonUnit. Then for every prime r | N:
- (i) the reduction mod r of every intermediate value equals the corresponding value of the
  same program run over F_r;
- (ii) every intermediate polynomial has the same degree over Z/N and over F_r.

Hence the computed pair reduces mod r to the correct Mumford representative of the
corresponding element of J_F(F_r).

*Proof.* The program is a sequence of ring operations, divisions with remainder by
polynomials whose top coefficient was checked to be a unit, inversions of such units, and
branches. Every branch (Euclid's "while r ≠ 0", Cantor's "while deg u > 2", the bits of r)
depends only on degrees and zero tests. The top coefficient of a normalised polynomial is a
unit mod N, hence nonzero mod r, so degrees agree. A coefficient ≡ 0 mod N is ≡ 0 mod r. Ring
operations, unit inversion, and division with remainder by a polynomial with unit leading
coefficient all commute with reduction mod r. Quotient and remainder are unique, and over Z/N
the quotient is just retained. Induction on the steps proves (i) and (ii).
- Over F_r, F is a squarefree monic quintic: F_5 = 1 and disc F is a unit by step 4. So
  Cantor's algorithm computes the group law of J_F(F_r) with identity (1, 0) (Cantor, Math.
  Comp. 48 (1987); Koblitz, J. Cryptology 1 (1989)).
- Over a field, no NonUnit can arise.
- Experimental check: G1b, 200 runs at 30-bit p, q, 0 mismatches (§8). ∎

**Proposition A (separation; genus-2 Lenstra 1987 Prop. (2.6); proved).** Let p ≠ q be primes
dividing N, and let τ satisfy:
- (i) p ≤ v;
- (ii) f_5 ≢ 0 and disc f ≢ 0 mod p;
- (iii) #J_f(F_p) is w-smooth;
- (iv) f_5 ≢ 0 and disc f ≢ 0 mod q;
- (v) ℓ ∤ #J_f(F_q), where ℓ := P^+(ord D_p) and D_p := D mod p ∈ J(F_p).

Then Algorithm C on τ outputs a proper divisor of N.

*Proof.*
1. By (ii) and (iv), the gcds in steps 2 and 4 are ≠ N. If either is ≠ 1, it is output.
   Otherwise step 5 runs.
2. D_p = [P − ∞] with P = (A, B) affine. It is nonzero, since a genus ≥ 1 curve has no
   degree-1 function, and deg u = 1. So ω := ord D_p > 1.
3. ω divides #J(F_p), which is w-smooth, so ℓ ≤ w.
4. For each prime r, r^{v_r(ω)} ≤ #J(F_p) ≤ (√p + 1)^4 < (⌊√v⌋ + 2)^4, so v_r(ω) ≤ e(r).
5. Let m := v_ℓ(ω) ≥ 1 and k_0 := ∏_{r<ℓ} r^{e(r)} · ℓ^{m−1}. Then ω ∤ k_0 and ω | k_0 ℓ.
   The schedule is increasing and repeats ℓ e(ℓ) ≥ m times, so it computes [k_0]D and then,
   at the next step, [k_0 ℓ]D.
6. Suppose no NonUnit arises up to and including [k_0 ℓ]D. By Lemma 3:
   - [k_0 ℓ]D reduces mod p to [k_0 ℓ]D_p = 0. Its u-part has degree 0 mod p, hence over Z/N.
     It is monic, so u = 1, and so [k_0 ℓ]D_q = 0 as well.
   - By (v), ℓ ∤ #J(F_q), so ord D_q is prime to ℓ and divides k_0 ℓ. Hence ord D_q | k_0 and
     [k_0]D_q = 0.
   - But [k_0]D_p ≠ 0 has u-degree ≥ 1, while [k_0]D_q has u-degree 0. This contradicts
     Lemma 3(ii) for the prefix ending at [k_0]D.
7. So a NonUnit arises, and by Lemma 2 its gcd is a proper divisor. ∎

Remark. Lenstra's (v) is the same condition, stated for #E(F_q). The point of (v) is that it
is a condition at q only. It concerns the full group order, so no information about ord D_q
is needed. With (v) available, uniformity of D_q in J(F_q) is not needed, and neither is any
upper bound on smooth Jacobian orders. That upper bound is not available in genus 2.

**Lemma 4 (supply at the second prime; proved from LPP II Thm 1.1).** Let q > 8100 be prime and
ℓ any prime. Let V_ℓ(q) be the set of τ ∈ F_q^7 such that (iv) holds and ℓ ∤ #J_f(F_q). Then
#V_ℓ(q) ≥ q^7 / (144,000 λ(q)).

*Proof.*
1. I_q contains ≥ ⌊q^{3/2}⌋ integers. At least (1 − 1/ℓ)q^{3/2} − 2 ≥ q^{3/2}/2 − 2 of them are
   prime to ℓ.
2. Removing LPP's ≤ q exceptions leaves ≥ q^{3/2}/3 values z, since q > 8100.
3. Each such z has ≥ q^{9/2}/(24,000 λ(q)) quintics f with nonzero discriminant and
   #J_f = z. These sets are disjoint for distinct z, which gives ≥ q^6/(72,000 λ(q))
   quintics f.
4. Each f comes from A(f) ≥ q − 4√q ≥ q/2 tuples (Lemma 1(c)). ∎

**Proposition B (per-curve success probability; proved from Lemma 1, Prop A, LPP II).** Let p
be the least prime factor of N, with p > 8100 and p ≤ v, and let q ≠ p be another prime
factor. Put S_w(p) := #{z ∈ I_p ∩ Z : P^+(z) ≤ w}. Then

  Pr[Algorithm C outputs a proper divisor] ≥ (S_w(p) − p) · (1 − 4p^{−1/2}) / (24,000 λ(p) p^{3/2}) · 1/(144,000 λ(q)).

*Proof.*
1. Let T_p be the set of residues t ∈ F_p^7 satisfying (ii) and (iii), and let
   E := {τ : τ mod p ∈ T_p and τ mod q ∈ V_{ℓ(τ mod p)}(q)}. Here ℓ(t) = P^+(ord D_p) is a
   function of t.
2. By Prop A, success ⊇ E.
3. By Lemma 1(b), Pr[E] = Σ_{t∈T_p} p^{−7} · #V_{ℓ(t)}(q) q^{−7}
   ≥ (#T_p/p^7) · min_ℓ #V_ℓ(q)/q^7.
4. Every non-exceptional w-smooth z ∈ I_p contributes ≥ p^{9/2}/(24,000 λ(p)) quintics. Each
   quintic contributes A(f) ≥ p − 4√p tuples. So #T_p ≥ (S_w(p) − p) p^{9/2} (p − 4√p)/(24,000 λ(p)).
5. Apply Lemma 4. ∎

**Cost of Algorithm C (proved).**
- Steps 1–4 take (log N)^{O(1)} bit operations. The discriminant is a 9×9 Sylvester
  determinant, computed with integer Bareiss.
- Step 5 performs ≤ Σ_{r≤w} e(r)(2⌈log₂ r⌉) ≤ 4π(w) log₂((⌊√v⌋+2)^4) + O(π(w)) group
  operations, which is O(π(w) log v).
- Each group operation is O(1) operations on polynomials of degree ≤ 10 over Z/N, plus O(1)
  extended gcds of degree ≤ 4 and O(1) inversions: O((log N)²) bit operations classically.
- Total per curve: O(w (log v)(log N)²/log w).
- The prime list up to w is sieved once per stage in O(w log log w).

### 3.2 Algorithm H (unknown p; standard doubling)

0. If N is even, output 2. Trial-divide by all primes ≤ 8100, and output any factor found.
   (N is assumed composite and not a prime power. Perfect powers can be detected and split in
   polynomial time by integer root extraction. A primality test can screen out prime N.)
1. For j = 14, 15, … (so v_j = 2^j > 8100):
   - set v_j := 2^j, w_j := ⌈L_{v_j}[1/2, 1]⌉, and h_j := ⌈L_{v_j}[1/2, 1 + ε(v_j)] (log N)³⌉,
     where ε(v) := (log log v)^{−1/2};
   - run Algorithm C h_j times with (v_j, w_j);
   - output the first divisor found.

## 4. Theorem R and proof

**Lemma 5 (smooth supply in I_p under RH; proved conditional on RH via Younis Thm 1.3 and CEP).**
Assume RH. Let w = w(p) satisfy log w = (s + o(1))√(log p log log p) for fixed s > 0. Then
S_w(p) ≥ p^{3/2} L_p[1/2, −1/s + o(1)] as p → ∞. In particular, S_w(p) − p ≥ ½S_w(p) for
large p.

*Proof.*
1. Let x := p² − ½p^{3/2} and h := p^{3/2}. Then (x, x + h] ⊂ I_p, so
   S_w(p) ≥ Ψ(x + h, w) − Ψ(x, w).
2. Since x < p², we have x^{3/4} < h ≤ x, so θ = 3/4 applies. For large p,
   (log x)^K ≤ w ≤ x^{1/2−ε} with K = K(3/4).
3. Younis Thm 1.3, or its α-free form (7.2), gives Ψ(x+h, w) − Ψ(x, w) ≥ ½ h Ψ(x, w)/x for
   large p. The error O(log(u+1)/log w) tends to 0, and α → 1.
4. Here u = log x / log w = (2/s + o(1))(log p / log log p)^{1/2} and
   log u = (½ + o(1)) log log p, so u log u = (1/s + o(1))√(log p log log p).
5. CEP then gives Ψ(x, w)/x = L_p[1/2, −1/s + o(1)].
6. Finally, p = p^{3/2} · p^{−1/2} = o(S_w(p)), because L_p[1/2, ·] = p^{o(1)}. ∎

**Theorem R (proved, conditional on RH and on Younis Thm 1.3 as stated).** Assume RH. On every
composite N that is not a prime power, Algorithm H outputs a proper divisor of N. Its expected
running time is

  L_p[1/2, 2 + O((log log p)^{−1/2})] · (log N)^{5+o(1)}   (p = least prime factor of N).

*Proof.*
1. **Correctness.** By Lemma 2, every output is a proper divisor.
2. Assume p > 8100; otherwise step 0 finds a factor. N is not a prime power, so it has a
   second prime factor q > p > 8100.
3. **The first adequate stage.** Let j* := min{j : v_j ≥ p}, so p ≤ v_{j*} < 2p. For j ≥ j*:
   - (i) of Prop A holds;
   - w_j ≥ w_{j*} = L_p[1/2, 1 + o(1)], and S_w(p) is increasing in w;
   - so Prop B and Lemma 5 (with s = 1) give a per-curve success probability
     π_j ≥ π* := L_p[1/2, −1 − δ(p)] · λ(q)^{−1}/C, with δ(p) → 0.
4. **Size of δ.** δ(p) collects:
   - the CEP o(1), which is O(log log u/log u) = O(logloglog p/ loglog p);
   - the change of variables in Lemma 5, which is of the same order;
   - λ(p), the constants, and the factor ½, which are L_p[1/2, O(log log p/√(log p log log p))].
   So δ(p) = O(logloglog p / log log p) = o(ε(p)). Since ε decreases and v_{j*} < 2p,
   ε(v_{j*}) ≥ ε(2p) ≥ 2δ(p) for large p.
5. **Each stage j ≥ j* succeeds with probability ≥ 1 − 1/e.**
   - (1 + ε(v))√(log v loglog v) is increasing in v, so h_j is nondecreasing in j.
   - λ(q) ≤ (log N)^{2+o(1)} and v_{j*} ≥ p, so h_j π* ≥ h_{j*} π*
     ≥ L_p[1/2, ε(2p) − δ(p)] (log N)^{1−o(1)}/C ≥ 1 for large N.
   - The h_j curves use independent randomness, so a stage fails with probability
     ≤ (1 − π*)^{h_j} ≤ e^{−1}.
6. **Cost per stage.** cost_j ≤ h_j · O(w_j log v_j (log N)²)
   = L_{v_j}[1/2, 2 + ε(v_j) + o(1)] (log N)^5.
7. **Growth of the cost.** cost_{j+1}/cost_j ≤ exp(3(√(log 2v loglog 2v) − √(log v loglog v)))
   ≤ exp(O(√(loglog v/log v))) ≤ e^{1/2} for large v.
8. **Summing.** E[time] ≤ Σ_{j<j*} cost_j + Σ_{m≥0} e^{−m} cost_{j*+m}
   ≤ (j* + Σ_m e^{−m/2}) cost_{j*} = O(log p) · L_{2p}[1/2, 2 + ε(p) + o(1)] (log N)^5.
   This is the claim, since L_{2p}[1/2, c] = L_p[1/2, c + o(1)]. ∎

**Corollary R1 (L_N form; proved).**
- If p ≤ N^β, then L_p[1/2, c] ≤ L_N[1/2, c√β], because log p ≤ β log N and
  log log p ≤ log log N.
- (log N)^5 = L_N[1/2, o(1)].
- If p stays bounded, the time is (log N)^{O(1)}.
- Hence the time is L_N[1/2, 2√β + o(1)] uniformly over N with p ≤ N^β.
- Lenstra–Pomerance (J. AMS 5 (1992), rigorous and unconditional) runs in L_N[1/2, 1 + o(1)].
  Dovetailing step by step costs a factor 2, and E[min] ≤ min E. So
  L_N[1/2, min(1, 2√β) + o(1)].

**Corollary R2 (complete factorisation; proved).**
- Recursively split composite pieces, testing each piece for primality (AKS, deterministic
  polynomial time) and for being a perfect power.
- There are ≤ log₂ N splits. Every composite divisor m of N with ≥ 2 distinct prime factors
  has least prime factor ≤ P_2(N), the second-largest distinct prime factor of N.
- Total expected time: L_{P_2(N)}[1/2, 2 + o(1)] (log N)^{O(1)}.

## 5. Audit table: every loss

| Step | Loss factor | Size | Source |
|---|---|---|---|
| trial division ≤ 8100, parity, perfect power | O(1)·(log N)^{O(1)} | negligible | — |
| tuple → curve distribution (A(f) weights) | 1 − 4p^{−1/2} | constant | Lemma 1, Weil |
| LPP II fibre lower bound at p | 1/(24,000 λ(p)) | L_p[1/2, o(1)] | LPP II Thm 1.1 |
| LPP exceptional set (≤ p values of z) | S_w(p) − p ≥ S_w(p)/2 | constant | Lemma 5 |
| LPP interval I_p (length p^{3/2}, versus 8p^{3/2} for Weil) | uses only part of the range | constant | — |
| smooth density in I_p, Younis/CEP | ½ · L_p[1/2, −1/s + o(1)] | main term | Younis Thm 1.3, RH |
| separation at q, condition (v) | 1/(144,000 λ(q)) | (log N)^{−2−o(1)} | Lemma 4 = LPP II at q |
| exponent length and group operations per curve | O(π(w) log v) | L_p[1/2, s + o(1)] | §3.1 |
| one Cantor operation over Z/N | O((log N)²) | polylog | — |
| unknown p: doubling v_j = 2^j, v_{j*} < 2p | O(log p) stages; L_{2p} vs L_p | L_p[1/2, o(1)] | §4 |
| extraction | any NonUnit is a proper divisor; Prop A forces one | none | Lemmas 2–3 |

No step costs a power of p. The only polylog-in-N term that is not L_p[1/2, o(1)] is
λ(q) ≤ (log N)^{2+o(1)}, together with the arithmetic cost. Both sit in the explicit
(log N)^{O(1)} factor, exactly as in LPP I Thm 1.1 and in Lenstra 1987.

Interval fit check (step 2 of the brief). LPP needs smooth numbers in I_p, which has length
p^{3/2} = x^{3/4+o(1)} around x ≈ p². Younis covers h ≥ x^θ for every fixed θ > 1/2, down to
y ≥ (log x)^{K(θ)}. Our w = L_p[1/2, 1] = exp((1/√2 + o(1))√(log x log log x)) is far above
(log x)^K. Both inputs fit with room to spare.

## 6. Optimising the constant

### 6.1 Smoothness cutoff
- With w = L_p[1/2, s]: per-curve work L_p[1/2, s], success probability L_p[1/2, −1/s].
- The total exponent s + 1/s is minimised at s = 1, giving 2.

### 6.2 Other genera, other groups (analysis model)
- **Genus g.** Orders are near p^g, in an interval of length ≍ p^{g−1/2} = x^{1−1/(2g)}.
  The smoothness probability at w = L_p[1/2, s] is L_p[1/2, −g/(2s)]. The cost s + g/(2s) has
  minimum √(2g), increasing in g: g = 1 gives √2, g = 2 gives 2, g = 3 gives √6 ≈ 2.449.
  - For g ≥ 3 the interval exponent 1 − 1/(2g) ≥ 5/6 is inside Younis's range. So the only
    missing input would be an LPP-type fibre lower bound, and even with it the constant is
    worse.
  - Higher genus is the wrong direction. A longer interval relative to the size does not help,
    because the cost depends only on log(order)/log p = g.
- **Genus 1 (ECM).** The constant would be √2, but the Hasse interval has θ = 1/2 exactly.
  - Younis needs θ > 1/2, and the RH saddle-point/contour method has the same barrier: the
    Perron error needs H ≥ x/h, and the horizontal and remainder terms need H ≤ x^{1/2−δ}.
  - Soundararajan (RH, 2010) reaches h ≍ u√x/ρ(u/2), but only with a lower bound of z x^{−ε}
    and only for y ≥ exp(5√(log x log log x)). That is not a density statement.
  - So ECM stays conditional on the classical short-interval conjecture (Lenstra 1987,
    Conj. (2.10)).
- **Split Jacobians / E_1 × E_2.** The order is #E_1·#E_2, a product of two Hasse-interval
  numbers. The smoothness probability is again L_p[1/2, −1/s], with no gain. Each factor
  still needs θ = 1/2 rigour.
- **Stage 2 / one large prime ≤ w².** The semismooth-to-smooth density ratio is (u log u)^{O(1)}
  = L^{o(1)}, so the constant does not change.
- **Order of D instead of #J.** This might be smaller, but no rigorous lower-bound handle
  exists. It is an L^{o(1)} effect heuristically.
- Conclusion: 2 is the floor of this route (proved for the analysis model: cost s + g/(2s)
  with g ≥ 2 forced by θ > 1/2). Heuristically 2 is also the true running time of HECM, since
  Younis gives an asymptotic, not only a lower bound.

### 6.3 Unconditional variant (V4)
LPP I Thm 1.3 (OCR of the Leiden scan; exponents partly garbled):
- c₁ = (1,980,000)^{1/3} (exponent as read);
- for x ≥ x₃, z = L_x[2/3, c₁] and x^{3/4} ≤ y ≤ x, there are at least
  y · exp(−(log x)^{1/3}(log log x)^{−1/3}) z-smooth integers in [x, x + y].

Exponent caveat: any lower bound y · exp(−o((log x)^{2/3})) suffices.
- Substituting this for Lemma 5 in Prop B and §4, with w = z at x = p², gives expected time
  L_p[2/3, 2^{2/3} c₁ + o(1)] (log N)^{O(1)}.
- This is LPP I Thm 1.1, whose algorithmic half (part III) was never published.
- Younis's unconditional Thm 1.1 (θ > 17/30, y ≥ exp(C(log x)^{2/3}(log log x)^{4/3})) gives
  exp(O((log p)^{2/3}(log log p)^{4/3})). That is slightly weaker in the log log power.

## 7. Literature and novelty

Read directly:
- **LPP I** (Phil. Trans. R. Soc. A 345 (1993) 397–408).
  - Thm 1.1 claims L_p[2/3, c₀](log n)² rigorously.
  - p. 398 states the conditional L_p[1/2, 2 + o(1)](log n)² under "a similar hypothesis"
    on short-interval smooth numbers, and the conjectural ECM bound L_p[1/2, √2].
  - So the constant 2 is LPP's.
- **LPP II** (Proc. LMS 84 (2002)): Thms 1.1–1.2. The analysis is deferred to a "forthcoming
  third paper".
- **Pomerance, "Elliptic curves: applications and problems"** (IMA Abel Conference, Jan 2011),
  slide 8: part III is "to appear (actually, 'to be written')". The method "stands (or will
  stand when completed)".
- **Lenstra 1987** (Ann. Math. 126):
  - Props (2.6)–(2.8) give the separation argument we transcribe. We read them via OCR of the
    wstein.org scan.
  - Conjecture (2.10) gives the conditional L_p[1/2, √2].
  - Lenstra notes that under GRH Prop (2.7) improves by a log log p factor. That is an
    L^{o(1)} effect.
- **Younis 2024**: Thms 1.1–1.3. The paper names no application to factoring.
- **Soundararajan 2010** (arXiv:1009.1591; RH): existence in [x, x + C(ε)√x] for y = x^ε,
  and the lower bound z x^{−ε} above. It says explicitly that this is "still not strong enough
  to be applicable to the analysis of Lenstra's algorithm".
- **Granville 2008 survey** ("Smooth numbers: computational number theory and beyond"), §4a:
  pre-2024 short-interval results.
  - Friedlander–Granville 1993: asymptotic for y ≥ exp((log x)^{5/6+o(1)}),
    √x y² exp((log x)^{1/6}) ≤ z ≤ x.
  - Harman 1991 and LPP 1993: lower bounds for y ≥ exp((log x)^{2/3+o(1)}).
  - Xuan 1999 (RH): existence in [x, x + √x (log x)^{1+o(1)}] for y = x^ε.
  - Hafner 1993 (RH, unpublished, an "almost all intervals" result per Soundararajan).
  - None gives, before 2024, an RH density lower bound at y = L[1/2] for intervals x^θ with
    θ < 1. This is why the corollary was not available earlier: Younis's RH theorem appears to
    be the first such input.

Searched and found nothing (see the literature-agent summary appended in §9):
- an RH- or GRH-conditional rigorous L_p[1/2, c] factoring bound;
- any citation of Younis applied to ECM, HECM or factoring;
- any published LPP part III.

Unconditional, almost-all:
- Rigorous "ECM works for most inputs" results exist (Pomerance 1987; Lenstra–Pomerance 1992;
  Hafner–McCurley 1989). They concern averages over auxiliary numbers or over n, not an
  L_p[1/2] bound for a given p.
- Matomäki–Radziwiłł (Ann. Math. 183 (2016)): almost-all short intervals for fixed u only.
- Jain (arXiv:2502.10530, per codex r18): y ≥ exp((log X)^{2/3+ε}).
- The obstruction at y = L[1/2] is the saddle point: 1 − α ≈ log u / log y ≫ (log x)^{−2/3},
  well outside the Vinogradov–Korobov zero-free region. Zero-density averaging over x does not
  obviously repair this, because the issue is ζ(s, y) near s = α for every x.
- So no unconditional (even almost-all-p) L_p[1/2] statement is known to us. Status of an
  almost-all-p theorem for genus 2: open; the r18 formulation of the missing lemma stands.

## 8. Numeric sanity check (pre-registered G1; post-hoc G1b)

Commands:
- `.venv/bin/python -m src.genus2rh.experiment_g1` (G1, 9 s)
- `.venv/bin/python -c "from src.genus2rh.experiment_g1 import main_b; main_b()"` (G1b, 30 s)
- Unit tests: `tests/test_genus2rh.py`, 6 pass. They cover group-law axioms; #J from exact
  F_p/F_{p²} point counts annihilating D and lying in the Weil interval; the monic model; the
  schedule; the discriminant; and proper divisors.

**G1** (40 pairs p < q in [101, 251], 30 curves each, v = 251, w = 50; 1200 runs):
- 358 runs COVERED. All 358 returned a proper divisor: **P1, 0 violations**.
- *But this is uninformative*: all 1200 runs returned a divisor.
  - 12 came at the f_5 stage and 18 at the discriminant stage (ACCIDENTAL. Pr[f_5 ≡ 0] = 1/r and Pr[disc ≡ 0] ≈ 1/r at each prime r,
    so the null rate is ≈ 2(1/p + 1/q) ≈ 2.4% per run; observed 30/1200 = 2.5%).
  - 1170 came in the main loop. For 786 of these, neither ord D_p nor ord D_q divides k.
    These are accidental coincidences in intermediate gcds. At p ≈ 200 there are hundreds of
    leading-coefficient tests per run, each hitting with probability ≈ 1/p.
- **P2 was never exercised**, because no run returned None. This is disclosed. Toy sizes cannot
  separate structural from accidental hits here.

**G1b** (post-hoc, disclosed; 20 pairs of 30-bit primes, 10 curves each, v = 2^30, w = 50):
- 200/200 runs returned None, as expected: 50-smooth Jacobian orders near 2^60 are rare.
- For all 200, the Z/N result reduced mod p and mod q equals the separate F_p and F_q
  computations: **P2, 0 mismatches**.

These checks confirm the implementation of Lemmas 2–3 and the Prop A implication. They say
nothing about the analytic inputs, which are asymptotic or need p > 8100.

## 9. Review log
(filled in after the GPT-6.1-sol referee pass `codex/genus2rh_1` and the Opus critic)

## Candidate card

- **ID / mechanism:** genus2rh-HECM. Genus-2 HECM (LPP) plus Younis's RH short-interval
  smooth density.
- **Target:** all composite non-prime-power N; time measured by the least prime factor p.
- **Representation:** random quintic y² = f(x) over Z/N with a known point; its Jacobian in
  Mumford form.
- **Extraction:** a nonunit leading coefficient in Cantor's algorithm, which occurs when
  [k]D vanishes mod p but not mod q.
- **Known lemmas:** LPP II Thm 1.1; Younis Thm 1.3 (RH); CEP; Weil; Cantor.
- **New lemma:** none analytic. The algorithmic accounting (Props A, B; Lemmas 3, 4) was
  unpublished (LPP III).
- **Cost:** L_p[1/2, 2 + o(1)] (log N)^{5+o(1)} (RH).
- **Closest prior:** LPP I's conditional remark, from which this differs only by naming the
  hypothesis (RH via Younis) and writing out the accounting; Lenstra ECM (conjectural √2).
- **Refutation criteria:** an error in Younis Thm 1.3; a misreading of LPP II's quantifiers
  (checked: ≤ p exceptions, every non-exceptional z, all quintics).
