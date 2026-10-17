# Phase 1 — mathematical findings (17 October 2026)

No polynomial-time factorizer is claimed. The report is the four-section Markdown response. Exact diagnostics are in phase1_experiments.py and phase1_results.json. Run: python3 -B phase1_experiments.py.

## Main correction

Formal univariate arithmetic-circuit PIT over Z/NZ is already randomized polynomial time, without knowing the factors: use random monic quotients of degree O(log degree), repeating polynomially often. Exact prior art: Arvind–Mukhopadhyay–Srinivasan, arXiv:0801.0514v1, Section 6, Algorithm 2, Lemmas 6.3–6.5, Theorem 6.6, crediting Agrawal–Biswas 2003. Thus binary support search needs a polynomial bound on pathwise ordinary-circuit section construction, not an additional unknown randomized PIT primitive.

Proof specialization: choose e>=2 with 2^e>=4D. If P is nonzero modulo N, choose p^a|N and the minimum coefficient valuation v<a; locally P=p^v Q with Q nonzero modulo p. Of p^e monic degree-e polynomials, at least p^e/(2e) are irreducible; at most D/e of these divide Q. Monic reduction commutes with multiplication by p^v and reduction modulo p. Hence a random quotient detects nonzero with probability at least 1/(4e). Unknown p appears only in the proof.

## Binary sections: exact linear-state obstruction

For s_t=[X^r]((1+X)^t mod (X^M-c)),
sum_t s_t z^t = z^r(1-z)^(M-1-r) / ((1-z)^M-cz^M).
The numerator and denominator are coprime when c is nonzero. For binary M=2^d, c=2, denominator degree is M over every odd residue characteristic. Thus every constant-coefficient linear realization generating all s_t has dimension at least M. This does not exclude nonlinear circuits or a special method for the single time t=N.

## Quotient ranks and a repaired false lead

For monic h of degree r, multiplication by f in F_p[X]/h has rank r-deg gcd(f,h). The coefficient pairing (u,v)->[X^(r-1)]uv mod h is nondegenerate, so the corresponding Hankel/Krylov matrix has the same rank. Full unit-pivot elimination over Z/NZ must encounter a useful nonunit if the two CRT ranks differ. Invertibility alone does NOT exclude nonunit entries or minors.

For odd N coprime to 3, W=X(X+1)(X^2+X+1) divides F_N over Z. For h=X^2-c, Res(h,Phi_3)=c^2+c+1. The toy norm splits N=35,c=2; N=77,c=2; N=143,c=3 are already found by this cheap resultant. We therefore classify Res(h,W) tests as baseline and inspect F only when W is a unit in the quotient.

After this repair, N=47959,r=2,c=46575 yields f=9209+16167X; coefficients and Res(h,W) are units, but det multiplication matrix =37810 mod47959 and gcd=199. CRT ranks are1 and2. This survives the cyclotomic baseline but is only a small success.

For N35, F35 mod5=2^5[X(X+1)Phi3^2]^5 and F35 mod7=5^7 W^7. Therefore avoiding W eliminates every CRT rank gap for ANY monic h, regardless of degree. Entries can still yield accidental factors.

For N>=3 the only cyclotomic factors of F_N are Phi2 (odd N, multiplicity1) and Phi3 (N=1 mod6, multiplicity2; N=5 mod6, multiplicity1). To exclude m>=4, evaluate at zeta=exp(2pi i/m): |1+zeta|^N>2 whereas |1+zeta^N|<=2. This is proved scaffolding, with no novelty claim.

## Discriminant grammar

P_K(T)=prod_(j=0)^(K-1)(T-j) has D_K=prod_(d=1)^(K-1)d^(2(K-d))=prod_(j=1)^(K-1)(j!)^2. Thus D_K is a unit modulo N iff least prime factor p>=K. D_(K+1)D_(K-1)=K^2 D_K^2. Cancellation after the first nonunit is invalid. The O(log K) shift-product description is a parameterized grammar, not an ordinary arithmetic SLP. BGS provides a sqrt(K) factorial-unit baseline; it does not automatically compute this superfactorial residue at that cost.

## Divided derivative obstruction

H_N=(1+X)^(N-1)-X^(N-1). For N=pq odd distinct primes and d=gcd(p-1,q-1), H_N has d-1 scalar roots in each component, so uniform scalar splitting probability is (d-1)(p+q-2(d-1))/(pq). The public root -1/2 is simple. Write N-1=2m, A=(1+X)^2,B=X^2. Then H_N=(2X+1) Q_N where Q_N=sum_(j=0)^(m-1) A^(m-1-j)B^j has a divide-free O(log N)-gate circuit. Q_N has exactly d-2 roots in each component. If d=2, Q_N(a) is a unit for every a modulo N.

## CRT masks

For odd ell coprime N, S_ell={s:s^2-4N is a square modulo ell, including zero}; |S_ell|=(ell+Legendre(N,ell))/2. Across L=prod ell the survivor density is exactly product |S_ell|/ell. A translated length-W interval has average W times that density, but the particular planted factor-sum interval has no proved bound. Efficient exact interval counting and few false survivors are separate missing lemmas. Conversion s->p has derivative magnitude p/(q-p), so Coppersmith transfer needs actual sufficiently narrow factor intervals and gap control.

## Independent norm-one construction

For quadratic algebra Z/NZ[sqrt(D)], let eps_p=(D/p),eps_q=(D/q), J=Jacobi(D/N), A=p-eps_p,B=q-eps_q,d=gcd(A,B). Sampling alpha uniformly among units and z=alpha/bar(alpha) gives uniform local norm-one groups. With m=N-J, gcd(m,A)=gcd(m,B)=d, and identity-in-exactly-one-component probability is d/A+d/B-2d^2/(AB). N35,D2 has probability0; N77,D2 has probability1/2. Varying D changes only four order pairs, so a generic success theorem needs more than twist selection.

## Scope of counters and tests

Primes used for the quotient scaling diagnostic are independently sampled by rejection at 8,12,16 factor bits; equal primes rejected; seed20261017. Candidate parameters use independent seed43019 and receive N,r,c only. Reference factors occur only in harness labels and field ranks. There are12 instances per size and48 parameter probes per instance,1728 total.

Counters instrument explicit quotient multiplication, GCD calls, and modular inversions; built-in modular exponentiation and GCD/inverse internal bit work are NOT expanded into ring operations. All raw records include N and diagnostic truth. No timing comparison or factoring speedup is claimed. The scalar/root/discriminant/CRT diagnostics deliberately enumerate small sets. The full seed harness explicitly materializes coefficients through degree1000; the candidate quotient degree is at most8, and the PIT demonstration uses degree at most12.

Separate mathematical critic and complexity auditor reviewed identities, state size, field assumptions, and the repair; final synthesis remains the main agent's responsibility.

Exact quadratic rank enumeration: after c(1-c)(c^2+c+1) is a unit, both local multiplication ranks are2 for EVERY c when N=101*113,107*193,1009*1481. For N199*241 the exact rank-gap probability is2/65. Therefore no positive uniform success bound exists for the repaired degree-two rank-only model. This does not exclude useful equal-rank nonunits or higher-degree quotients.
