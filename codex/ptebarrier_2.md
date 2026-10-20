# ptebarrier round 2 (ADVERSARIAL REFEREE): try to break the rigidity barrier

You are an adversarial referee (analytic number theory / Diophantine approximation / computational
number theory). Your job is to find errors, gaps, overclaims and counterexamples. Do not be polite.

Read `notes/ptebarrier_results.md` in this repository (and, if useful, `src/ptebarrier/rigidity.py`,
`notes/systeq_results.md` §5–6, `notes/ptecons_results.md` §3). The repository is read-only for you;
put scratch work (Python is fine; `.venv/bin/python` has sympy, python-flint, gmpy2) in the current directory.

Check, in this order, and give a verdict for each (CORRECT / CORRECT WITH FIXES / WRONG / UNCLEAR):

1. **Lemma R (§5.1)**, every step. In particular:
   (a) the Mahler-measure height bound on G̃ = G/g at integer points |a| ≤ B;
   (b) the "at most L far roots" step and the treatment of multiplicities;
   (c) the degree/leading-coefficient step when deg Q ≠ deg Q';
   (d) the Newton/Vandermonde step with complex roots of Q, Q' (is J computed right? is μ ≠ 0?);
   (e) the case analysis (good/bad roots), especially Case 2, and whether the final inequality
       M ≥ K + 2L + 4j + 3 really covers all cases.
   Try hard to construct counterexamples just inside the claimed regime, e.g. by computer search for
   pairs (Q P, Q' P') with small low-degree difference, or by structured constructions (Chebyshev/
   Dickson compositions, PTE chains, Prouhet–Thue–Morse, products with common factors).
2. **Theorem W (§5.2)**: the reduction to Lemma R (height h of G = Q_{A'} r_A − Q_A r_{A'}, degree
   K = k + j₂, multiplier degree j = j₁ + j₂), and the conclusion "same rational function ⇒ same
   relation" in both views (residue μ, and integer evaluation m). Any hidden assumption (units,
   zero divisors, r_A(m) = 0, content of polynomials)?
3. **Corollary P (§5.3)** and **Prop S (§5.4)** (Lagrange interpolation with |F(m)| ≥ N; the
   "reduced relation" cancellation; the PTE consequence M − k ≥ d + 2ρ − 1).
4. **Prop NX / ONE-RING (§5.5)**.
5. **§3 (Vinogradov) and §4 (determinant method, k = 1)**: is the k = 1 argument correct (linear
   components, components not over Q, multiplicities), and is the Binyamini–Cluckers–Novikov
   citation/statement right (absolute constant for plane curves? exponent of log?)?
6. **§6 (Theorem PB, the L-range translation, and "what it rules out")**: is the claimed scope
   honest? Is the model 𝓜_sys (exact identities, one presentation, polynomially small small side) a
   faithful formalisation of systeq §6's open escape, or does it quietly exclude the natural
   versions? Are O1–O4 the right list of what remains open, and is anything listed as closed actually
   open (or vice versa)? Is the View-A claim in §6.2 (|r(m)| ≥ |m|^{deg r}/2 ⇒ regime) right?
7. **Prior art**: is Lemma R (or Theorem W / Prop S) already known (PTE constant lower bounds,
   Rees–Smyth, Caley, Borwein; "polynomials close to split polynomials"; Mason–Stothers-type
   statements; Lenstra lacunary results)?

Label every finding FATAL / MAJOR / MINOR with a concrete fix. Give an overall verdict
(ACCEPT / ACCEPT WITH FIXES / MAJOR REVISION / REJECT). Keep it under ~3500 words.
