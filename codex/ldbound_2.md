# ldbound_2: adversarial referee round on Theorem LD* for design (a'''_k)

You are a hostile, meticulous referee in analytic number theory. The repository
/Volumes/git/factoring-connection is read-only for you (scratch: your cwd). Your job is to BREAK the
claimed proof, not to be encouraging. Rigour over speed.

## Read
1. /Volumes/git/factoring-connection/notes/ldbound_results.md — the claims under review, in full.
2. /Volumes/git/factoring-connection/notes/pcest_results.md — the proved inputs it relies on
   (Lemma GB, Lemma 1.2, Lemma 1.3, Theorem FM and its Steps 1–4, Lemma 8.5', Lemma Q, Lemma DS,
   Theorem R, Lemma Kern, Prop 4.2, Theorem 4.3, Theorem A). These passed two reviews; you may still
   flag anything the new claims use incorrectly.
3. Primary sources as needed: Montgomery–Vaughan, "The large sieve", Mathematika 20 (1973) (Brun–
   Titchmarsh with constant 2) and the additive large sieve (K − 1 + Q^2); Lee–Venkatesan arXiv 1805.08873;
   Harper arXiv 1208.5992.

## The claim in one paragraph
Design (a'''_k) adds a uniform random prime w in (W, 2W], W = N^{s+1/(4k)+eta0/2}, and takes
F(g0, g1) = N w (Lemma GB for N w), with height H = N^h, h >= s + 1/(4k) + eta0. Then the forced residue
F^(P) mod t (t = |G(P)|) is c_P w with c_P independent of w, so by Brun–Titchmarsh its law is <= C/phi(t)
pointwise (Lemma W). Kernel + AM–GM + enlargement + averaging over w reduce E[N(N-1)] to an all-residue
pair count Q(t l q) = #{(F, F') smooth in [-X*, X*]: F = F' mod t l q}, summed over the rational fibres
(Lemma DS bounds the partner's rational smoothness along each fibre) — Steps 0–7 of §5. The all-residue
count is bounded at precision L^{o(1)} (Theorem BDH_fam, §4) by additive Parseval plus: Lemma DS for
conductors <= B^{C0}; the additive large sieve for larger conductors, using the divisibility sparsity
nu(m') <= 2 X'_+ q / m_1 (class B) and the average over the special-q primes l (class C). Conclusion:
PC <= L^{o(1)}, so Theorem A holds for (a'''_k) conditional on LV6 (+ HECM, irreducibility transfer),
with theta''' = 3/(d+1) + d/(4k), same infimum (1 + 3/(d+1))/sqrt 2.

## What to do
For EACH of: Lemma GB_w, the sizes paragraph, Lemma 1.2_w, coprimality (1.1), the FM re-check (§2),
Lemma W, Corollary W, the Parseval identity, Theorem BDH_fam (A), (B), (C) and (4.1), Theorem LD*
Steps 0–7, Theorem A''', Props 7.1–7.3 and the constants table, give a verdict CORRECT / CORRECT WITH
FIXES / WRONG / UNCLEAR with precise reasons and line references, and the repair if one exists.

Attack in particular:
- Hidden dependence on w anywhere other than r0 (windows, filter, candidate set, counts, sizes, Kern).
  Is E_w legitimately taken inside the pair sum? Is anything conditioned on w?
- Step 2's AM–GM and Step 3's enlargement: do they lose more than a constant at the main-term scale?
- Step 5: is the fibre parametrization right (direction (a, b), step -Delta in G(P'))? Is Lemma DS
  applicable exactly as used (interval inside [X B^{-C0}, X], length >= X B^{-C0}, modulus <= B^{C0})?
- BDH_fam: check every inequality numerically in exponents of N and L; the large sieve's (K + Q^2) term,
  the multiplicity nu(m'), whether each m' is counted with the right multiplicity, and whether the l-average
  is legitimately used given that other factors (#P_l, U_l, V_l) depend on l. Are there conductors or
  moduli that fall between the classes (A), (B), (C)? Non-reduced classes? F = 0? Negative F?
- The final power count V^4/Y^2 ~ S^4, and the comparison with FM's lower bound (E N)^2.
- FM for (a'''_k): is Step 2(d)'s window argument really fine with |F^| up to M w V^d? Does Lemma 8.5' still
  apply with the stated ranges? Does anything in FM use F(g0, g1) = N exactly?
- Does F(g0, g1) = N w break anything downstream: Theorem K / field selection (hypD), the relation map,
  LV §6 (square roots, characters, the g1 column), irreducibility, the "every nonunit is output" rule?
- Is the claimed constant table right (theta''' = 3s + d/(4k) at h = s + 1/(4k); min k)?
- Is there a cheaper design (smaller h) achieving the same, or is Prop 7.3's barrier right?
- Literature: is Theorem BDH_fam (or something implying it) already published? Is Lemma W's use standard?

End with: an overall verdict on "Theorem LD* for (a'''_k) is proved", the list of required fixes, and
whether Theorem A''' is a correct conditional reduction on LV6 (+ the listed caveats).
