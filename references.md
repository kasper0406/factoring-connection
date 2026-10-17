# References

Primary sources consulted in this project, with the version consulted. Provenance codes:
**[A: agent]** the named agent read the primary text; **[G]** checked by GPT-6.1-sol only (prior-art map or codex rounds),
not re-read by an agent; **[G, partial]** GPT-6.1-sol reports partial access; **[cited]** cited but not checked by anyone.
Section pointers name the note that records what was checked.

## Brief reading list (BRIEF §8)

- **[AKS]** M. Agrawal, N. Kayal, N. Saxena, *PRIMES is in P*. Annals of Mathematics 160 (2004) 781–793. Journal PDF. Lemma 2.1, §2. [G] `prior_art_map.md` §1.
- **[BGS]** A. Bostan, P. Gaudry, É. Schost, *Linear Recurrences with Polynomial Coefficients and Application to Integer Factorization and Cartier–Manin Operator*. SIAM J. Comput. 36(6) (2007) 1777–1806. Author-hosted final PDF. Lemmas 3, 4, 12, 13; Theorems 5, 8, 11, 15. [A: factC] `trackC_tal_audit.md`; [G] `prior_art_map.md` §2.
- **[Harvey]** D. Harvey, *An exponent one-fifth algorithm for deterministic integer factorisation*. arXiv:2010.05450v1 (12 Oct 2020). Thm 1.1, Prop. 2.5, Lemma 2.4. [G] `prior_art_map.md` §3, `codex_r1_tal_bgs.md`.
- **[HH]** D. Harvey, M. Hittmeir, *A log-log speedup for exponent one-fifth deterministic integer factorisation*. arXiv:2105.11105v1 (24 May 2021). Thm 1.1. [G] `prior_art_map.md` §4.
- **[Sparse]** J. van der Hoeven, G. Lecerf, *Sparse polynomial interpolation: Faster strategies over finite fields*. Final version, 2 May 2024 (texmacs.org). Prop. 2, Alg. 4, §3.2. [G] `prior_art_map.md` §5.
- **[Tal]** Y. Tal, *Computing Modular Factorials Below the Square-Root Barrier*. ECCC TR26-211, posted 25–26 Sep 2026, 69-page PDF. §§1–5, 6.1–6.2, 8, App. A, H. [A: factC] `trackC_tal_audit.md`; [G] `codex_r1_tal_bgs.md`.
- **[Coppersmith]** D. Coppersmith, *Small Solutions to Polynomial Equations, and Low Exponent RSA Vulnerabilities*. J. Cryptology 10 (1997) 233–260. Journal PDF, §§10–11, Thms 2, 4, 5. [G] `prior_art_map.md` §7. Thresholds re-derived via [May10] by trackE.
- **[Regev]** O. Regev, *An Efficient Quantum Factoring Algorithm*. arXiv:2308.06572v3 (7 Jan 2024). Thm 1.1, Lemma 4.4. [G] `prior_art_map.md` §8.
- **[Pilatte]** C. Pilatte, *Unconditional correctness of recent quantum algorithms for factoring and computing discrete logarithms*. arXiv:2404.16450v3 (24 Aug 2026); Forum Math. Pi 14 (2026) e5. Thms 1.1, 3.18. [G] `prior_art_map.md` §9.
- **[EG]** M. Ekerå, J. Gärtner, *Extending Regev's factoring algorithm to compute discrete logarithms*. arXiv:2311.05545v2 (30 Jan 2024). Lemma 2. [G] `prior_art_map.md` §10.
- **[Gao]** Y. Gao, Y. Feng, H. Hu, Y. Pan, *On Factoring and Power Divisor Problems via Rank-3 Lattices and the Second Vector*. arXiv:2512.19076v1 (22 Dec 2025). Thm 1.1, Lemma 4.4. [G] `prior_art_map.md` §11 (statements only; proof not reproduced).
- **[Cleve]** R. Cleve, *The query complexity of order-finding*. arXiv:quant-ph/9911124v1 (30 Nov 1999). Thms 1–2. [G] `prior_art_map.md` §12.

## Lattice decoding and partial information

- **[May10]** A. May, *Using LLL-Reduction for Solving RSA and Factorization Problems: A Survey*. In *The LLL Algorithm* (2010) 315–348. Author's 34-page preprint (archived `crypto.rub.de` copy), §§2–5, Thms 1–3, §5.1. [A: trackE] `trackE_analysis.md` header.
- N. Howgrave-Graham (lemma as stated in [May10] Thm 2). [A: trackE, via May10].
- N. Heninger, H. Shacham, *Reconstructing RSA private keys from random key bits*. CRYPTO 2009. [cited] exact theorem not re-checked, `trackE_analysis.md` §4.6.
- M. Herrmann, A. May, *Solving linear equations modulo divisors*. Asiacrypt 2008. [cited] bound derived, not re-verified, `trackE_analysis.md` §3.1.
- O. Goldreich, D. Ron, M. Sudan, *Chinese remaindering with errors* (2000); D. Boneh, *Finding smooth integers in short intervals using CRT decoding*, STOC 2000. [cited] not tested, `trackE_analysis.md` §3.3.
- D. Boneh, G. Durfee, N. Howgrave-Graham, *Factoring N = p^r q for Large r*. CRYPTO '99, LNCS 1666, 326–337. Author upload. Thm 2, Lemma 1. [G] `prior_art_map.md` §35.
- A. May, M. Ritzenhofen, *Implicit Factoring: On Polynomial Time Factoring Given Only an Implicit Hint*. PKC 2009. Thms 5, 7. [G] `prior_art_map.md` §37.
- F. Sica, *Factoring with Hints*. arXiv:0912.1585v6 (rev. 8 Aug 2017; PDF dated 26 Oct 2018). Thm 1. [G] `prior_art_map.md` §38.

## Straight-line programs, factorials and many-root polynomials

- R. J. Lipton, *Straight-line complexity and integer factorization*. ANTS I, LNCS 877 (1994) 71–79. Two-page publisher preview only; full theorem UNVERIFIED. [G, partial] `prior_art_map.md` §13.
- R. J. Lipton, *Factoring and Factorials*, blog post, 23 Feb 2009. Thm 1. [G] `prior_art_map.md` §14.
- M. Shub, S. Smale, *On the intractability of Hilbert's Nullstellensatz and an algebraic version of "NP ≠ P?"*. Duke Math. J. 81 (1995) 47–54; manuscript dated 21 Jun 1994. Main Theorem, Thm 2. [G] `prior_art_map.md` §15; statement via Bürgisser [A: factC].
- P. Bürgisser, *On Defining Integers And Proving Arithmetic Circuit Lower Bounds*. Computational Complexity 18 (2009) 81–103; ECCC TR06-113 (2006). Thm 1.1 (ECCC) / Main Thm 1.2 (journal). [A: factC, ECCC version] `trackC_analysis.md` §3; [G] `prior_art_map.md` §16.
- Q. Cheng, *On the Ultimate Complexity of Factorials*. TCS 326 (2004) 419–429. Author manuscript. Thm 1 (conditional). [G] `prior_art_map.md` §17.
- Q. Cheng, *Straight Line Programs and Torsion Points on Elliptic Curves*. Computational Complexity 12 (2003) 150–161. Author manuscript. [G] `prior_art_map.md` §18.
- P. Dutta, G. Jindal, A. Pandey, A. Sinhababu, *Arithmetic Circuit Complexity of Division and Truncation*. CCC 2021, LIPIcs 200, Art. 25. Thms 6/26, 30, 32. [A: factC] `trackC_analysis.md` §2 M3; [G] `prior_art_map.md` §24.
- E. Allender, P. Bürgisser, J. Kjeldgaard-Pedersen, P. B. Miltersen, *On the Complexity of Numerical Analysis*. SIAM J. Comput. 38(5) (2009) 1987–2006. Author manuscript. Props 2.3, 2.4. [G] `prior_art_map.md` §19.
- P. Koiran, S. Perifel, *The Complexity of two Problems on Arithmetic Circuits*. TCS 389 (2007) 172–181; manuscript 24 Aug 2007. Numbered positive-characteristic theorem UNVERIFIED. [G, partial] `prior_art_map.md` §20.
- H. Fournier, G. Malod, S. Mengel, *Monomials in arithmetic circuits: Complete problems in the counting hierarchy*. STACS 2012; manuscript 11 Apr 2012. Thms 3.1–3.3. [G] `prior_art_map.md` §21.
- N. Kayal, C. Saha, *On the Sum of Square Roots of Polynomials and Related Problems*. ECCC TR10-189 (2010). Thms 1.5/4.5, 2.1. [G] `prior_art_map.md` §22.
- D. Perrucci, J. Sabia, *Real roots of univariate polynomials and straight line programs*. J. Discrete Algorithms 5 (2007) 471–478. Thm 7. [G] `prior_art_map.md` §23.
- C. Umans, S. Wang, *A number-theoretic conjecture implying faster algorithms for polynomial factorization and integer factorization*. arXiv:2511.10851v1 (13 Nov 2025). Thms 5.2, 5.5. [G] `prior_art_map.md` §32.

## Identity testing, quotient rings, generic models

- M. Agrawal, S. Biswas, *Primality and Identity Testing via Chinese Remaindering*. JACM 50(4) (2003) 429–443; manuscript 21 Feb 2003. Thm 4.8, Lemma 4.7. [G] `prior_art_map.md` §28.
- V. Arvind, P. Mukhopadhyay, S. Srinivasan, arXiv:0801.0514v1, §6 (Algorithm 2, Lemmas 6.3–6.5, Thm 6.6). [G] `codex_r1_mechanisms.md` A1; `codex/PHASE1_NOTES.md`.
- D. J. Bernstein, *An exposition of the Agrawal–Kayal–Saxena primality-proving theorem*, 10 Aug 2002. Not fetched; UNVERIFIED. [G, partial] `prior_art_map.md` §29.
- D. Aggarwal, U. Maurer, *Breaking RSA Generically is Equivalent to Factoring*. EUROCRYPT 2009. Thm 1, Lemmas 6–8. [G] `codex_r2_generic_ring.md` §3.
- T. Jager, J. Schwenk, *On the Analysis of Cryptographic Assumptions in the Generic Ring Model*. IACR eprint 2009/621, full version 25 Jan 2012; J. Cryptology (2013). Thms 2–3. [G] `codex_r2_generic_ring.md` §3, `codex_r4_sos_classgroup.md` §1.
- I. Damgård, M. Koprowski, *Generic Lower Bounds for Root Extraction and Signature Schemes in General Groups*. EUROCRYPT 2002. Thm 1, Cor. 1. [G] `codex_r2_generic_ring.md` §3.
- D. Boneh, R. Venkatesan, *Breaking RSA may not be equivalent to factoring*. EUROCRYPT 1998. Thms 3.1, 3.4. [G] `codex_r2_generic_ring.md` §3.
- U. Maurer, S. Wolf, *Lower bounds on generic algorithms in groups*. EUROCRYPT 1998. Thm 10. [G] `codex_r2_generic_ring.md` §3.
- D. R. L. Brown, *Lower Bounds for Straight Line Factoring*. IACR eprint 2010/047 (revised). §5. [G, with caveats] `codex_r2_generic_ring.md` §3.

## Number theory: orders, characters, sieves, class groups

- J. M. Pollard, *Theorems on factorization and primality testing*. Proc. Cambridge Philos. Soc. 76 (1974) 521–528. Publisher extract only. [G, partial] `prior_art_map.md` §36.
- H. W. Lenstra Jr., *Factoring integers with elliptic curves*. Annals of Math. 126 (1987) (description of Pollard `p−1`). [G] `codex_r3_critic.md` §2.
- H. C. Williams `p+1` method, as implemented in GMP-ECM `pp1.c`. [G] `codex_r1_mechanisms.md` F1.
- C. Pomerance, sieve account and CFRAC implementation study (math.dartmouth.edu/~carlp). [G] `codex_r3_critic.md` §2, `codex_r2_generic_ring.md` §4.
- D. W. Boyd, G. Martin, M. Thom, *Squarefree values of trinomial discriminants*. LMS J. Comput. Math. 18 (2015) 148–169, §5, Lemma 5.2. [A: repair] `trackBD_analysis.md` §2 B5.
- G. Martin, A. Parvardi, *Subproducts of small residue classes*. Canad. Math. Bull. 65(1) (2022) 1–8, Cor. 3 (cites Burthe). [A: repair] `trackC_analysis.md` §3 O4.
- J. L. Massey, *Shift-register synthesis and BCH decoding* (1969). [G] `codex_r3_critic.md` §3.
- Z. Kelley, arXiv:1602.00208, Thms 1.1, 2.3 (lacunary root bounds); Canetti–Friedlander–Konyagin–Larsen–Lieman–Shparlinski (root bounds). [G] `codex_r3_critic.md` §3.
- D. A. Burgess, character-sum bounds. [cited] UNVERIFIED citation, `trackC_analysis.md` §3 E8.
- M. Hittmeir, hyperbolic sieve, arXiv:2205.10074, §5. [G] `codex_r2_generic_ring.md` E2.
- M. Hittmeir, *Smooth Subsum Search*. arXiv:2301.10529v2 (11 Jul 2023). [G] `prior_art_map.md` §33.
- P. Stănică, E. Mulder, M. Hittmeir, *Multi-target hyperbolic sieves and elliptic trace obstructions*. arXiv:2606.13018v2 (5 Oct 2026). Thm 3.12, Cor. 3.16. [G] `prior_art_map.md` §34.
- Z. Oznovich, B. L. Volk, *On Deterministically Finding an Element of High Order Modulo a Composite*. arXiv:2506.07668v3 (11 Oct 2025). Thm 1.1. [G] `prior_art_map.md` §30.
- D. Harvey, M. Hittmeir, *Deterministic methods for finding elements of large multiplicative order*. arXiv:2601.11131v2 (5 Jun 2026). Thm 1.1. [G] `prior_art_map.md` §31.
- OEIS A020639, *Least prime factor of n* (S. Sykora contribution, 29 Jul 2017). [G] `prior_art_map.md` §26.
- L. M. Adleman, K. S. McCurley, *Open problems in number theoretic complexity, II* (1994), C11/O11a–b. [G] `codex_r4_sos_classgroup.md` §1.
- R. Schoof, *Elliptic Curves over Finite Fields and the Computation of Square Roots mod p*. Math. Comp. (1985), §4. [G] `codex_r4_sos_classgroup.md` §1.
- M. O. Rabin, J. O. Shallit, *Randomized algorithms in number theory*. CPAM (1986). Theorem numbering UNVERIFIED. [G, partial] `codex_r4_sos_classgroup.md` §1.
- E. Bach, J. Shallit, *Algorithmic Number Theory*, vol. I (1996), §§5.7–5.9, 7.1–7.2. [G, contents only] `codex_r4_sos_classgroup.md` §1.
- A. Rybalov, *On generic complexity of the problem of representation of natural numbers by sum of two squares* (2020). [G] `codex_r4_sos_classgroup.md` §1.
- P. Stevenhagen, arXiv:1806.06250 (Rédei matrices; Thms 2.1, 3.1, eq. (16)); corrected Rédei reciprocity, Math. Proc. Cambridge Philos. Soc. (doi 10.1017/S0305004121000335). [G] `codex_r4_sos_classgroup.md` §3.
- P. Milović, 4-ranks of `Cl(−4pq)` (UCL Discovery `dpq8rank.pdf`), Prop. 1. [G] `codex_r4_sos_classgroup.md` §3; consequence verified numerically by `src/repair/verify_r4.py`.
- C. P. Schnorr, H. W. Lenstra Jr., class-group factoring (1984); H. W. Lenstra Jr., C. Pomerance, JAMS 5(3) (1992). [G] `codex_r4_sos_classgroup.md` §4.
- Negative Pell equation `x² − Ny² = −1` (odd continued-fraction period of `√N` ⇒ solvable ⇒ all odd prime factors ≡ 1 mod 4); discussion cited by the GPT-6.1-sol r5 review: math.toronto.edu/~eknight/Pell.pdf. [G] `codex/r5_report_review.out` item 9; used as the positive control in NEXT_STEPS item 1.
- R. Schoof, *Quadratic fields and factorization* (reneschoof.github.io `mcpaper.pdf`), pp. 265–266; J. E. Gower, S. S. Wagstaff Jr., *Square form factorization*, §4.5. [G] `codex_r4_sos_classgroup.md` §4.

The full GPT-6.1-sol prior-art map, with the exact hypotheses and complexity variables of each entry, is
`notes/prior_art_map.md`. "VERIFIED" there means GPT-6.1-sol checked the statement against the source, not that the proof
was verified.
