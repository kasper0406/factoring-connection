You are a research mathematician (computational number theory). Be concrete, terse, label every claim
PROVED / KNOWN (cite) / CONJECTURED / REFUTED. Work in your scratch dir; you may write and run Python
(python3; install cypari2/sympy with pip --user if needed) to check claims on small inputs.

Setting: N = pq, p<q odd primes, unknown. Goal of the project: a poly(log N) (or sub-L[1/3]) classical
factoring algorithm. Already established (do not redo): genus characters of Cl(-4kN) = QR_N; square roots
in Cl(D) factor N; ambiguous-form routes need the group order or relations; a fixed discriminant gives O(1)
bits which can be brute-forced; class numbers to precision <1 need ~sqrt(N) terms of an approximate
functional equation.

Answer these four questions, each in <= 40 lines, with a verdict line first.

Q1 (Kronecker-Hurwitz mod 2^k). sum_{t^2<4N} H(4N-t^2) = 2 sigma(N) + (terms computable from N) gives
p+q, hence a factor, from class numbers of a FAMILY of ~sqrt(N) discriminants. Is there ANY identity
(Eichler-Selberg trace formula at other weights/levels, Cohen-Zagier, Hecke-operator traces mod 2^k,
Nicolas-Serre / Bellaiche-Nicolas mod-2 Hecke algebra of level 1, mod-2^k Hecke algebra) that gives
sigma(N) mod 2^k, or tr T_N mod 2^k on some space, from poly(log N) quantities each computable natively
from N? In particular: in the Nicolas-Serre framework T_p acting on F_2[Delta] depends on Frob_p in the
pro-2 extension of Q unramified outside 2. For N=pq, T_N = T_p T_q: is any coefficient of T_N (e.g. the
action of T_N on Delta^k mod 2 for explicit k) both (a) computable natively from N in poly(log N), and
(b) not a function of N mod 2^j for j=O(1)? Coefficient of q^N in Delta^k mod 2 counts representations by
diagonal forms sum 2^i m_i^2 - say exactly which counts, and whether any of them is known to be computable
in poly(log N) without factoring (I believe not: it is a lattice count mod 2).

Q2 (QR oracle => factoring?). An oracle for the Legendre symbol (a/p) for arbitrary a mod N (equivalently
a QR_N oracle). State the best known result: is any poly-time or sub-L[1/3] factoring algorithm known using
it? Also: oracle for (r/p) restricted to primes r <= (log N)^2 (determines p under GRH): any known
efficient decoder? Same for a quartic-residue oracle and for "QR in Z[sqrt a]/N on one split component"
(Redei symbols [a,b,p]). Cite (Goldwasser-Micali open problem, Jager-Schwenk generic model, Boneh-Venkatesan,
Damgard Legendre PRF, Seurin 2-Phi-hiding, Hittmeir).

Q3 (p mod 2^j from 2-power class groups). Tower results: for N=1 mod 8, 8-ranks of Cl(-4N), Cl(-8N) give
p mod 16 in one cell. Does any family of 2^j-ranks / Redei-type symbols of Cl(-dN), d public, determine
p mod 2^j for growing j? (Coppersmith needs j ~ (log2 N)/4.) My expectation: no -- p mod 2^j is Frobenius
in Q(zeta_{2^j}), and 2^j-rank layers instead carry spin / non-abelian Frobenius (Koymans-Milovic 16-rank
of Q(sqrt(-2p)), FIMR spin). Confirm or refute with literature and a small computation.

Q4 (native route of your own). Propose ONE concrete poly(log N) candidate for a single 2-power class-group
bit beyond genus of a discriminant attached to N (e.g. r4 Cl(-4N), h(-4N) mod 8) that does not obviously
need sqrt mod N, relations, group order or N^c work, and locate exactly where it fails. Test it on small N
against PARI (qfbclassno / quadclassunit) and report agreement rates on >= 200 N.

Write your final answer as markdown.
