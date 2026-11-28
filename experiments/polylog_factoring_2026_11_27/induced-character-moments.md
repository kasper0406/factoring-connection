# The first full-period moment that separates the hidden primes

Let $\chi$ be a primitive odd quadratic character of short conductor $f$, with $N=rd$ a product of distinct odd primes and $\chi(r)=1$, $\chi(d)=-1$. The sign condition $\chi(N)=-1$ is public. Put $m=fN$ and $\psi(a)=\chi(a){\bf1}_{(a,N)=1}$. For the raw moments $M_j=\sum_{a=1}^{m}\psi(a)a^j$, the first three vanish and the cubic moment is

$$
M_0=M_1=M_2=0,\qquad M_3=NC(1-r^2)(1+d^2),\qquad
C=\sum_{a=1}^{f}\chi(a)(a^3-f^2a).
$$

Thus a **single modular reading** of the full unit-period sum gives a factor on balanced inputs. The earlier reduction used a product on an interval of length $\sqrt{2N}$. [Earlier split-prime product](collective-norm-product.md) · [Executable check](check_induced_character_moments.py)

## Elementary proof and decoding

Put $u_j=\sum_{a=1}^f\chi(a)a^j$. Odd symmetry gives $u_0=0$ and $u_2=fu_1$. Summing powers on $h$ successive periods gives
$$
S_0(h)=0,\quad S_1(h)=h u_1,\quad S_2(h)=fh^2u_1,\quad
S_3(h)=h C+f^2h^3u_1.
$$
Möbius inversion imposes the unit condition as $M_j=\sum_{e\mid N}\mu(e)\chi(e)e^j S_j(N/e)$. Every term proportional to $\sum_{e\mid N}\mu(e)\chi(e)$ vanishes. The cubic term left over gives the displayed formula. Equivalently $C=fB_{3,\chi}=-3fL(-2,\chi)$. This constant is nonzero by the primitive functional equation and the Euler product at 3. Computing $C$ directly takes one period of the short character. [Generalized Bernoulli numbers](https://dlmf.nist.gov/24.16.E10) · [Functional equation](https://dlmf.nist.gov/25.15.E5) · [Special values of Hurwitz zeta](https://dlmf.nist.gov/25.11.E14)

Suppose $\max(r,d)<2\min(r,d)$. Reading $M_3\bmod |C|N^2$, dividing by $|C|N$, adjusting the sign of $C$ and subtracting one gives $h=d^2-r^2\bmod N$. Since $|h|<3N/2$ for the actual difference, at most three lifts need testing. For each, check whether $\sqrt{h^2+4N^2}$ is an integer and whether half the sum and difference are squares. A final division verifies the recovered factor. If $(C,N)=1$, $M_3\bmod N^2$ alone suffices using $C^{-1}\bmod N$. The lower digit modulo $N$ always vanishes. More generally the sum of $\psi(a)g(a\bmod N)$ over one period is zero for **every** function $g$ on the units modulo $N$, including reciprocals, by CRT and $\sum_{a\bmod f}\chi(a)=0$.

## Full-period analytic sums

Inducing the character gives
$$
\frac{L(s,\psi)}{L(s,\chi)}=(1-r^{-s})(1+d^{-s}).
$$
At $s=1$ the ratio is $1+(r-d-1)/N$, so absolute error below $1/(2N)$ suffices to recover $r-d$ by rounding. At $s=0$, differentiation gives $L'(0,\psi)=2L(0,\chi)\log r$. These identities use only the missing Euler factors and need no zero-location hypothesis. [Induction formula](https://dlmf.nist.gov/25.15.E4)

Odd symmetry and the partial fractions of cotangent give
$$
L(1,\psi)=\frac{\pi}{2m}\sum_{a=1}^{m-1}\psi(a)\cot(\pi a/m).
$$
Grouping by $a\bmod f$ leaves a sum over primitive $N$th roots within each class. With $\zeta=\exp(2\pi i/f)$ this becomes
$$
L(1,\psi)=-\frac{\pi i}{fN}\sum_{u\in(\mathbf Z/f)^\times}
 \chi(u)\frac{\zeta^u\Phi_N'(\zeta^u)}{\Phi_N(\zeta^u)} .
$$
Indeed the roots in one class are $\xi\omega$ with $\xi=\zeta^{aN^{-1}}$ fixed and $\omega$ ranging over primitive $N$th roots. The sum of $(\xi\omega+1)/(\xi\omega-1)$ is $\varphi(N)-2\xi^{-1}\Phi_N'(\xi^{-1})/\Phi_N(\xi^{-1})$. The constant term cancels, and $\chi(-N)=1$. At $f=4$, $r\equiv1$, $d\equiv3\bmod4$, one has $\Phi_N(i)=1$ and $i\Phi_N'(i)=\varphi(N)/2+i(N-1+r-d)/2$. Obtaining the factor-dependent derivative is the expensive step in this rearrangement. [Cotangent expansion](https://dlmf.nist.gov/4.22.E3) · [Earlier cyclotomic reductions](bounded-quadratic-roots.md)

The finite Fourier transform does not move the exceptional indices to shorter frequencies. With the positive-exponent convention, CRT gives $\widehat\psi(k)=\chi(N)\tau(\chi)\chi(k)c_N(k)$. Here $c_N(k)$ is a Ramanujan sum, equal to 1 when $(k,N)=1$. At every index below the smaller prime, the original sequence equals $\chi(a)$ and its Fourier transform equals the public background $\chi(N)\tau(\chi)\chi(k)$. Whenever $N\nmid k$ and the Fourier coefficient differs from this background, $\gcd(k,N)$ is a proper factor. The identity allows aggregate methods but does not give a short list of revealing coefficients. [Ramanujan, *On certain trigonometrical sums and their applications*, 1918](https://ramanujan.sirinudi.org/Volumes/published/ram21.pdf) · [Gauss sums](https://dlmf.nist.gov/25.15.E6)

At $N=77,f=8$, $C=72$. The exact cubic moment is $-33{,}264{,}000$. Its residue $33{,}264\bmod 72N^2$ recovers 11. The checks recover the positive-character prime in 63 balanced fixtures across six discriminants. Exact finite-field tests cover 1,677 Fourier coefficients and three cotangent identities. The direct summation uses $O(fN)$ terms and is a reference evaluator. Evaluating the unrestricted polynomial moments from one short period is inexpensive. Imposing coprimality by the displayed inclusion-exclusion formula requires the unknown factors. [Checks](check_induced_character_moments.py)
