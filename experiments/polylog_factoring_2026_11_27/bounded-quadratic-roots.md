# Recovering mixed roots from partial oracle answers

A univariate small-root lattice can finish the decoder once the oracle has confined a root to an interval of length $N^{1/2-\varepsilon}$ for any fixed $\varepsilon>0$. The same holds for an affine progression with that many terms and step coprime to $N$. Thus the nonadaptive quarter-window reduction needs only $(1/2+\varepsilon)\log_2 N+O_\varepsilon(1)$ scaled disk queries, plus its reference. The implemented five-dimensional lattice uses a wider safety margin and needs $(3/5)\log_2 N+O(1)$ scaled queries. [Quarter-window reduction](nonadaptive-quarter.md)

## Explicit lattice

For a monic quadratic $f(x)\bmod N$ and an integer bound $X$, use the coefficient vectors after substituting $x\mapsto Xx$ in
$$N^{m-i} f(x)^i,\qquad xN^{m-i}f(x)^i\quad (0\le i<m),\qquad f(x)^m.$$
Each polynomial is divisible by $N^m$ at a root of $f\bmod N$. The basis is triangular by degree and has determinant $N^{m(m+1)}X^{m(2m+1)}$. LLL and Cauchy–Schwarz give a nonzero combination whose absolute value for $|x|\le X$ is at most
$$\sqrt{2m+1}\,2^{m/2} N^{m(m+1)/(2m+1)}X^m.$$
When this is strictly less than $N^m$, every bounded modular root is an ordinary integer root of the combination. Factor that integer polynomial and verify each candidate against the original congruence. Fixed $m$ permits $X\le c_m N^{m/(2m+1)}$ for a positive constant $c_m$. The code uses $m=2$ and $X=\lfloor N^{2/5}/4\rfloor$, with direct tests for intervals of at most two integers.

For an interval or progression $u=r+ht$, substitute into $u^2-1$ and divide by the unit $h^2\bmod N$ to make the quadratic monic. In the quarter-window reduction, partial binary digits give an interval for the lower mixed root. Stop generating oracle queries once the known interval has at most $X+1$ integers. Decode by the lattice and take gcds of the roots minus 1. For the affine CRT decoder, the same argument permits stopping the auxiliary prime list once its product exceeds the window length divided by $X+1$. If the trial retained a singleton, the CRT residues then leave an arithmetic progression short enough for the lattice. Extra candidates can always be verified. [Affine CRT decoder](affine-parity.md)

## Why short windows alone do not complete the reduction

A uniform random translate retains each mixed root with probability $w/N$, even for a fixed unit multiplier. Thus with $w=O(\sqrt N)$ the probability of retaining any mixed root is at most $2^{\omega(N)}w/N=N^{-1/2+o(1)}$. Here $\omega(N)=O(\log N/\log\log N)$ follows from the product of the smallest distinct primes. This bounds the probability of successful isolation before the lattice is applied. The lattice gives a native polynomial-time evaluator for the contents of already narrow affine windows, while the reduction needs wide windows or a source of partial coordinates.

## Checks

The code matched exhaustive root enumeration on 300 bounded progressions, including a repeated-prime input. It also split three inputs end to end after querying the shortened quarter-window transcripts. The 19-bit input 471,311 took 13 scaled disk queries and one reference. The code uses SymPy for lattice reduction and integer polynomial roots. Run `python bounded_quadratic_roots.py`.
