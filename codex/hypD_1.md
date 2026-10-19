# hypD_1 — blind research round: rigorous relation supply for Dixon descent (L[1/2, c<1] factoring)

You are a research number theorist. Read-only repo at /Volumes/git/factoring-connection (scratch: cwd).
Context files (read): notes/lemmaM_results.md (Theorem DD, regimes (a),(a'), Hypothesis D, three missing
ingredients), notes/nfsfruit_results.md §1 (Lee–Venkatesan randomized NFS, arXiv 1805.08873).

GOAL: a rigorous randomized algorithm factoring EVERY N in expected time L_N[1/2, c] with c < 1
(beating Lenstra–Pomerance 1992), or the strongest provable partial result and the sharpest barrier.

Setting (Dixon descent): R = r^2 mod N uniform square; split R = u/v (mod N), |u|,|v| <= sqrt(N)/2; the
large primes l in (B, B1] of u, v are "descended" by special-q relations in an LV-randomized field
A = (m, f = fhat_{N,m} + (x - m)R(x)), coefficients of R uniform of height H. Fruitfulness is free
(Theorem DD). What is open is relation supply.

Please work on these questions, labelling every claim PROVED / KNOWN (with exact citation and exact
ranges of moduli, residues, smoothness y) / CONJECTURED / FALSE:

Q1 (literature, exact statements). Give the precise current statements, with exact uniformity ranges, of:
 (i) Harper's Bombieri–Vinogradov for smooth numbers (arXiv 1208.5992), including the shape of the error
     (is it Psi(x,y)(e^{-cu/log^2 u} + y^{-c}) + Q sqrt(Psi) log^{7/2} x?);
 (ii) Drappeau 2015 (Compositio, Fouvry–Iwaniec type for friable integers), level and residue conditions;
 (iii) Pascadi 2025 (arXiv 2505.00653) level 5/8 (or better), residue conditions;
 (iv) individual-modulus results: Granville 1993, Soundararajan 2008, Harper 2012 (q <= y^{4 sqrt e - delta}),
     exact y-ranges, and whether lower bounds Psi(x,y;q,a) >> Psi_q(x,y)/phi(q) are known for larger q;
 (v) BV for moduli restricted to SPARSE sets (e.g., smooth moduli, LV Lemma 8.5; multiples of a fixed l);
 (vi) smooth values of binary forms / polynomials at small smoothness (Balog–Blomer–Dartyge–Tenenbaum,
     Lachand, Harper "smooth values of polynomials", Dartyge) — exact smoothness ranges;
 (vii) 2018–2026: any rigorous L[1/2, c<1] factoring or DL in prime fields, any rigorous descent analysis
     (Gordon 1993, Schirokauer, rigorous index calculus, Enge–Gaudry, Kleinjung–Wesolowski), LV follow-ups.

Q2 (key idea to test). Claim: ingredients 2 and 3 of Hypothesis D collapse into ONE averaged statement.
 For a field A let G_A = set of primes l in a dyadic range I of large primes that are "descendable" in A.
 For any set S of integers in [1, X], X = sqrt(N)/2, Pr_R[R = u/v with u, v in S] >= (const)*2^{-omega} *
 (|S|/X)^4 (positivity: sum_{u,v in S} 1[u/v in (Z/N)^{x2}] = 2^{-omega} sum_chi |sum_{u in S} chi(u)|^2 >=
 2^{-omega}|S|^2, plus a gcd/multiplicity argument), and 2^omega = L^{o(1)} WLOG after Pollard–Strassen
 removes prime factors < exp((log N)^{1/2}(log log N)^{1/4}). Taking S_A = {m * l_1...l_k : m B-smooth,
 l_i in G_A}, the descent in field A needs only density g_A = |G_A|/pi(I) >= exp(-o(log L / k)).
 So Hypothesis D follows from: Pr_{A, l uniform in I}[l descendable in A] >= exp(-o(log L/k)).
 Check this rigorously (including multiplicities, coprimality, the square-class condition, the split law,
 and codex lemmaM_2's objection that two field types may descend complementary prime types).

Q3 (first moment). With l averaged, the rational-side moduli are t = a - bm in T = {l*s: l in I, s B-smooth}.
 Does LV's Lemma 8.5 (BV restricted to smooth moduli with saving rho(Q,y)) extend to moduli in T, with
 the numbers counted being B-smooth and the moduli B1-smooth (two different smoothness bounds)?
 What exact range of Q relative to Psi(x,B) does it need, and does LV's "goodness" (Psi(F,B;r,a) =
 (Psi_r/phi(r))^{1+o(1)}) suffice for an L^{o(1)}-loss first moment at L[1/2] scale?

Q4 (second moment / existence). From a first moment E_{A,l} N_l(A) we need Pr_{A,l}[N_l(A) > 0] — a
 second moment (pairs of special-q pairs in one field) with loss <= exp(o(log L/k)), k ~ (log N)^{1/4}.
 Is there ANY known way to bound Pr_A[F_A(a,b) and F_A(a',b') both B-smooth] (LV random R; joint image
 lattice of (<c,v>,<c,v'>) has index gcd of 2x2 minors) with such loss? Or a design avoiding second moments
 (e.g., many independent fields MNFS-style; choosing which primes to descend; one-large-prime partials)?
 Alternatively: is there a design where per-target success follows from FIRST moments only?

Q5 (designs and constants). For each design you consider, compute the best L[1/2, c] constant attainable
 using only PROVED inputs, and identify the single missing estimate. Include: (a') coefficient inflation;
 QS-style descent (x^2 - N = l*c); Dixon-only descent with Vallee's 2/3 sampler; a "spanning" design where
 cheap non-fruitful LV relations plus ONE fruitful Vallee/Dixon relation suffice (needs a rank statement);
 L[1/2] designs where log-scale losses are affordable.

Be concrete and rigorous; give formulas. Do not pad. Final answer <= 4000 words.
