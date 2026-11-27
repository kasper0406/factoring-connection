# modrec: can a lattice-count residue mod 2^k evade the hull barrier? (2026-11-27)

**Verdict: no poly(log N) evaluator was found. Both the main hypothesis and my own idea are killed for
uniform mechanisms.** Barvinok mod 2 is also killed in the binomial-denominator model (codex). Pre-registration: `experiments/modrec/prereg.md`. Code: `experiments/modrec/`.

## Yield ledger

| Mechanism | What it computes exactly | Cost / growth | Status |
|---|---|---|---|
| Hull (continued-fraction Stern–Brocot) descent, Pick bookkeeping mod 2^k (`hull_count.py`) | #{(x,y) in [0,X]^2 : F(x,y) <= T} for definite F, from N alone. The floor8 short-square query S(N) mod 8 is 1 + sum_{j<=8} (C(jN) - C(jN-1)) | Agreed with the factor-backed `lattice_sum` on 184/184 cases (two fixtures, plus 30–54 bits, D = -23 and -4). Hull vertices fit N^{0.333} for both forms. At 54 bits: 1.06M vertices and 18.6 s (D = -23); 35-bit fixture 15k vertices and 0.16 s | Exact **N^{1/3}** factor-blind evaluator for the open floor8 query (experimentally supported). Reducing the accumulator mod 2^k leaves the visited nodes unchanged, because the control flow depends on exact geometry (proved) |
| Mod-2 pruning by involutions of descent nodes | A node region bounded by an arc of {F <= T} is mapped to itself by an affine map only if that map fixes the centre and is an automorph of F. Aut(x²+xy+6y²) = {±1}, which never preserves a first-quadrant cap. For x²+y² only the D4 reflections qualify, so just the caps straddling an axis or the diagonal | 0 certifiable prunes for D = -23. O(1) per quadrant for D = -4 | Refuted for this model (proved). Symmetries of a cap's lattice set alone cannot be certified without building its hull |
| Uniform mod-2 reciprocity for incomplete quadratic floor sums | Let F(p) = sum_{y<p/2} floor(2y²/p). For primes p = 3 mod 4, p > 3, h(-p) = 3 - 2F(p) mod 4. Equivalently, ((p-1)/2)! = ±1 mod p (Mordell 1961). Proof: floor(2y²/p) is odd iff y² mod p > p/2, then apply Dirichlet's class number formula | Identity held on all 9,005 primes below 2·10^5; Mordell checked on 1,135 primes below 2·10^4 (`classno_parity.py`) | Proved and verified. A uniform mod-2 descent for incomplete quadratic counts at the prescribed half-period cutoff would compute h(-p) mod 4 for primes in poly(log p). We know of no such algorithm, even with no factoring involved. This does not cover a method restricted to rank-2 (ellipse) instances |
| Whole-disk count mod 16 | C(T) = 1 + 4(Sq(T) + 2G(T)) mod 16. G counts n <= T with r_2(n)/4 = 2 mod 4, i.e. n = p^e·(square or twice a square) with p = 1 mod 4, e = 1 mod 4 | Verified for every T <= 10^6 (`disk_mod16.py`) | Proved and verified. The mod-16 disk count is a split-prime counting parity: the Polymath4/Tao–Croot–Helfgott class, where π(x) mod 2 is known only in x^{1/2-c}. Mod 8 is trivial; mod 16 is already prime-counting. That hardness is an analogy, not a reduction |
| **Own idea: modular forms mod 2^k** | The D = -23 hidden bit is a(N) mod 2 for η(z)η(23z), with a(pq) = a(p)a(q) from the splitting of x³-x+1. Euler's pentagonal theorem gives a(N) = #{u² + 23v² = 24N, u = v = -1 mod 6} mod 2 | Matched the q-expansion for all N < 3000. Matched the cubic-splitting formula on 17,590 semiprimes below 10^5 (`eta23_mod2.py`) | Killed. The mod-2 reduction maps the representation bit at N to a representation parity of the same type at index 24N: no shrinkage. The cutoff 1/(1-q) has every coefficient odd. Serre's density theorem says almost all coefficients of an integral-weight form vanish mod 2^k, so the cutoff cannot be made modular mod 2^k. Only the analytic approximate functional equation remains, and it is not exact mod 2^k |
| Zagier/windmill-type involutions | On {x² + 4yz = N}, swap fixes the representations x² + 4y² = N. The windmill map fixes factorisations N = y(y+4z) | Classical | Parity only: both fixed sets are even for N = pq = 1 mod 4. Information mod 4 needs the whole reduction cycle, i.e. the class group of Q(√N), which costs L[1/2] |
| Codex modrec_1 (blind, session 01a21e40-dafa-75a0-8803-390c0bb6758c) | Exact halving descent for H(T,X) = sum_{d<=X} floor(T/d) mod 2: H(2n+b,X) = H(n,X/2) xor (one threshold test per odd d <= X). I re-checked it on 3,000 random cases. Its F_2 rational GF obstruction: odd primes in (X/2, X] survive in the boundary's reduced denominator, so any sum of binomial-denominator fractions needs Omega(X) exponent bits | X - log X threshold tests in total, i.e. Θ(√T) at X = √T, which is worse than the hull. The GF lower bound covers only the Barvinok-type binomial-denominator model | Killed (proved). It also confirmed that the octant increment at N = pq is 2 or 0, so mod 4 is needed, not parity. Code: `experiments/modrec/codex_modrec_1/` |

## Reformulation worth keeping
For x²+y² the short-square query is S(N) = 1 + 4·[N is a sum of two squares] mod 8 (floor8 already
records this). So a fast evaluator would decide quadratic residuosity of -1 mod N, i.e. recognise Blum
integers. For D = -23 it decides a hidden-class bit (a(N) mod 2 above). The cutoff count matters only
through differences at jN and jN-1. Region methods must therefore resolve the level sets
F = jN, which are representation numbers.

## Sharpest open step
A descent that works only for definite rank-2 forms, at a single prescribed cutoff, and that uses
the congruence N | F. None of the barriers above covers it: the parabola barrier needs rank 1, and
the prime-count analogy needs arbitrary T. Such a method would decide QR(-1 mod N), a classical
open problem. Astra and I found no candidate mechanism.
