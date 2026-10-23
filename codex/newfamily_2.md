# newfamily_2: ADVERSARIAL referee round

You are a hostile but fair referee (number theory, complexity, statistics). The repo
/Volumes/git/factoring-connection is read-only; write scratch files only in the current directory. Put your FULL
report in your final message (not only in a scratch file).

Read: `notes/newfamily_results.md` (the paper under review), `notes/newfamily_prereg.md`, the code
`src/newfamily/{recip,modsym,battery,harness,diag}.py`, `tests/test_newfamily.py`, and your own earlier blind report
`codex/newfamily_1_report.md`. Context: `notes/factorbits_results.md` §3.2 (Theorem PT, PT-m),
`notes/selfred_results.md` (Theorem L), `notes/nonsmooth_results.md` §5.4.

Attack, in this order, and for each give a verdict CORRECT / FIXABLE (state the fix) / WRONG (counterexample):
1. **Theorem RA** (§2.1): the proof (conductor, m = generator of 𝔣₀ ∩ Z, positivity at real places, N coprime to m
   but maybe not to the discriminant of K — does Art need N·O_K coprime to anything else?), the claim that it covers
   de Boer–Pagano power-residue symbols at a rational modulus, the claim "L/Q need not be abelian", and
   **Corollary RA-PT** (pure trees + symbol nodes ⊂ PT-m). Is the corollary's model sensible (e.g. symbol of c_v where
   c_v mod N is computed but c_v itself is a huge N-independent integer)? Prior literature for RA?
2. **Prop DS/DS2/DSF** (§2.3): the decomposition, the BHK formula with (−1)^r term, the congruence claims (12Ns mod 3,
   mod 8, mod N), the factoring proof of DSF (k_p congruence at p; 2aX ≡ (p+q)(a²+1) mod N; edge cases p = 3, a² ≡ −1,
   p+q ≥ N?), and the formula 6N·U ≡ (a²+1)φ(N)/(2a) (mod N). Is DSF known (e.g. in the literature on Dedekind sums
   and factoring, or "Dedekind sums of composite moduli")?
3. **Prop SAL** (§2.4): the K(a,b;N²) closed form for composite N (twisted multiplicativity + prime-square evaluation
   — any missing Jacobi/ε factor for p ≡ 3 mod 4?), the decoder's correctness and the heuristic failure analysis,
   and whether "constant precision" is honestly supported. Is the overall Salié/Kloosterman-at-square ⇒ factoring
   statement known?
4. **Class E** (§2.5) and the "random argument" reduction: is the conditional statement right? Is the equidistribution
   input (Bettin–Drappeau etc.) available uniformly in the factorisation of an individual modulus?
5. **Battery** (§3): validity rules, degenerate statistics, the n = 32 size artefact explanation (Chebyshev-type
   finite-prime-set imbalance), the nonunit accounting (golden-ratio argument), MDE, and whether the conclusion is
   stated no more strongly than warranted.
6. **Candidate table** (§1): any verdict you dispute; any important candidate missing (be wild: you may add up to 5
   new candidates with a verdict).
7. **Overall**: is the verdict ("no fourth family; class E is the only computable non-aggregate bin; most promising
   open candidate = class E at structured arguments") justified? What single claim is the weakest?

Output: numbered findings with severity (MAJOR / MINOR / NIT), each with a concrete fix, then a one-paragraph overall
recommendation (ACCEPT / MINOR REVISION / MAJOR REVISION).
