# What theta, Hecke and finite Poisson transforms preserve

The fixed-cap disk can be written exactly as a coefficient problem for a marked theta series. This does not give a faster evaluator. At the low coefficients needed by the decoder, the ordinary Hecke formula returns the original coefficient at index $Nj$. Finite Poisson transformation moves the four square-root slopes to four dual slopes, including the two hidden ones. These identities reject two direct modular-inversion shortcuts while identifying the extra structure a new transform must use.

## Four index-$N$ lattices

Let

$$Q_0(x,y)=x^2+bxy+dy^2,\qquad 4d-b^2=\Delta,$$

and let $c^2\equiv-\Delta\pmod N$. For each of the four roots $r$ of $r^2\equiv-\Delta\pmod N$, put

$$s_r=(r-b)/2\pmod N,
\qquad
F_r(u,v)=N u^2+(2s_r+b)uv+\frac{Q_0(s_r,1)}N v^2.$$

The last coefficient is integral and $F_r$ has discriminant $-\Delta$. The map

$$(u,v)\longmapsto(Nu+s_rv,v)$$

is a bijection from $F_r(u,v)=j$ to the root line inside $Q_0(x,y)=Nj$. For $1\le j\le A<p$, every solution has $y$ invertible modulo $N$. If $y$ vanished at either prime, the quadratic congruence would also force $x$ to vanish there, contradicting the valuation of $Nj$. The four root lines therefore partition all representations.

The two public roots $r=\pm c$ give forms equivalent to $Q_0$. The other two give inverse forms with the same ordinary theta series. Calling either hidden form $F$, the exact coefficient identity is

$$R_{Q_0}(Nj)=2R_{Q_0}(j)+2R_F(j),\qquad 1\le j\le A.$$

The hidden form is the square of the prime ideal class. It is generally not the prime form $[p,c,q]$. For $N=28829=127\cdot227$ and $j=7$, substituting $[127,339,227]$ predicts four representations, while the exact left side and the correct hidden-square contribution are both zero.

## Why the Hecke index does not shrink

For the weight-one theta series with quadratic character $\chi_{-\Delta}$, the coefficient formula at an index $N$ coprime to the level is

$$[q^j](T_N f)=
\sum_{e\mid(j,N)}\chi_{-\Delta}(e)[q^{Nj/e^2}]f.$$

Here $j\le A<p$, so $(j,N)=1$. The sum has one term and becomes

$$[q^j](T_N f)=[q^{Nj}]f.$$

Thus writing the desired high coefficient as a low coefficient after $T_N$ changes its label and does not compute it. Building $T_N$ from $T_pT_q$ uses the unknown factorization. A direct double-coset or coefficient construction must still be charged. This is the same issue noted for large-index Hecke routines in the [weighted disk investigation](../ellipse_aggregate/README.md#aggregating-all-norm-indices).

The absolute-coordinate scalar contains still more information than the ordinary theta coefficient. On the $r$-line the original coordinate is

$$L_r(u,v)=2Nu+(2s_r+b)v.$$

The full signed disk sum of $|L_r|$ is four times the positive-quadrant scalar. The first derivative of an ordinary Jacobi theta series instead weights a vector by $L_r$. It vanishes coefficientwise because $(u,v)$ and $(-u,-v)$ cancel. The absolute value is the operation that prevents that cancellation. An ordinary holomorphic theta derivative therefore cannot replace the required marked absolute sum.

## The finite Poisson transform is self-dual in the hidden slopes

Let $\mathcal R=\{r:r^2=c^2\pmod N\}$ and, on the full residue square, use the line-multiplicity function

$$M(w,y)=\sum_{r\in\mathcal R}{\bf1}_{w=ry}.$$

Its unnormalised finite Fourier transform is

$$\widehat M(a,t)
=\sum_{r\in\mathcal R}\sum_{y\bmod N}e_N((ar+t)y)
=N\#\{r\in\mathcal R:ar+t=0\}.$$

For a unit $a$, this is nonzero precisely when $-t/a$ is one of the same four square roots. The two mixed slopes survive unchanged as dual slopes. On the fixed-cap disk every nonzero $y$ is a unit, so the line-multiplicity function agrees with the quadratic-congruence indicator there. Multiplication by the elliptical cutoff convolves this four-line spectrum with the transform of the boundary. It does not turn the hidden slopes into public data.

This is an exact finite identity. It is not a lower bound for every Poisson or automorphic algorithm. It shows that a useful transformation must also aggregate the sharp boundary or the absolute weight. Transforming only the complete congruence retains the original hidden-root problem.

## Checks and cost conclusion

The implementation checked 90 balanced branches, 474 levels and 1,280 source representations. It verified the four lattice bijections, the theta coefficient identity, the marked absolute-coordinate identity and the dual-slope formula. It also records the $N=28829$ prime-form counterexample.

```sh
python experiments/polylog_factoring_2026_11_27/theta_hecke_transform.py
```

All form construction and coefficient decoding use integers with $O(\log N)$ bits once a root is supplied. Supplying a mixed root is factoring-equivalent. Computing the high coefficients by direct enumeration, or evaluating every dual frequency, costs a power of $N$ in the current methods. No polynomial-bit evaluator follows. The remaining admissible target is a transform that handles the fixed ellipse boundary and the absolute coordinate jointly, without constructing the mixed slopes or expanding a power-sized frequency set.
