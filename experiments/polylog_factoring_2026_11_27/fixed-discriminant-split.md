# Factoring reduction at one discriminant and cutoff

The fixed-window reduction can use primitive forms with one common discriminant, one leading coefficient and one cutoff for each odd modulus. Each shifted query takes one signed count modulo 4. One reference query cancels the horizontal axis throughout, including for nonsquarefree moduli. The forms can vary in class. [Fixed-window bisection](fixed-window-split.md)

## Common parameters

For $N\ge5$ odd, put $M=(N-1)/2$ and $K=\lfloor M/2\rfloor$. Take a power of two $L>16N^2$ and set $C_0=L^2(K+1)^2$. Choose the even integer $D\in[C_0,C_0+2N)$ congruent to $-L^2$ modulo $N$. Fix
$$ B_0=\left\lceil\frac{C_0+(LK+2N)^2+2N}{N}\right\rceil.$$
For window start $a$, put $s=a+K+1$ and take the odd integer $j\in[0,2N)$ congruent to $Ls$ modulo $N$. Query the form
$$Q_j(x,y)=(Lx+jy)^2+D y^2.$$
Its leading coefficient is $L^2$, and its discriminant is $-4L^2D$. Its last coefficient $j^2+D$ is odd, so the form is primitive because $L$ is a power of two. The whole query description has $O(\log N)$ bits.

## Exact rows and the reference residue

Write $H=NB_0$. At $y=1$ all $|x|\le K$ satisfy $Q_j(x,1)\le H$. For $|x|\ge K+1$, we have
$$Q_j(x,1)\ge C_0+(L(K+1)-2N)^2>H.$$
The last inequality follows because the difference of the two squared bounds is $L(2K+1)(L-4N)>3N$. Also $H<2C_0$, since $L^2(2K+1)>4NLK+4N^2+3N$. Thus every row with $|y|\ge2$ lies outside the disk. The cutoff is below $\sqrt{4L^2D/3}$ because $B_0<2C_0/N<L^2(K+1)$. At $y=0$ the form is $L^2x^2$, invariant under every shift.

Modulo $N$, dividing the form by the unit $L^2$ gives $(x+sy)^2-y^2$. The rows $y=\pm1$ together contribute twice the count of roots in $[a+1,a+2K+1]$. Subtract twice the known-root count and, when needed, twice the extra endpoint indicator as in the fixed-window construction. The corrected signed residue is $R(a)=U+2W(a)\pmod4$, with a single unknown axis contribution $U$.

The midpoint $a_0=\lfloor M/2\rfloor$ has known mixed parity $W(a_0)=0$ unless the directly testable integer $a_0+1$ is a mixed root and $M$ is odd. That exceptional case splits immediately by gcd. Otherwise $R(a_0)$ supplies the common baseline. Comparing all subsequent residues with it permits the same $O(\log N)$ bisection. Primality and perfect-power testing handle the leaves. This establishes the reduction without a disk evaluator whose bit cost is polynomial.

## Checks

The independent signed-disk reference verified 6,478 windows on 158 odd moduli, including the common cutoff, common discriminant, primitive forms and invariant axis residue. The reduction split all 84 inputs with at least two distinct prime factors. Run `python fixed_discriminant_split.py`.
