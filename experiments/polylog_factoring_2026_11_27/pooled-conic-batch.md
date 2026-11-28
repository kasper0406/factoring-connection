# Sharing a norm-product polynomial across rows

The earlier conic search tests $\gcd(N,i^2+Aj^2)$ at $i\ge0$, $j\ge1$, $i^2+Aj^2\le H=\lfloor\sqrt{8N}\rfloor$, with $A=x^2+1$. All rows can share the same polynomial evaluations. This preserves the exact set of candidate norms and takes
$$
 N^{o(1)}\bigl(1+N^{1/4}/A^{1/4}\bigr)
$$
bit operations for $A\le H$. At $A$ comparable to $\sqrt N$ it retains the earlier $N^{1/8+o(1)}$ bound. Which inputs have a witness at a given $x$ is unchanged. [Earlier row-by-row algorithm and coverage bounds](large-conic-batch.md) · [Pooled implementation](pooled_conic_batch.py)

For a dyadic interval of $w$ columns shared by $b$ rows, take a power-of-two stride $h$ within a constant factor of $\max(1,\sqrt{w/b})$. For the row indices $j$ form
$$
 P(Y)=\prod_j\prod_{0\le k<h}\bigl((Y+k)^2+Aj^2\bigr)\pmod N.
$$
Evaluate at the $w/h$ block starts. Only a nonunit block needs individual gcds, at most $hb$ of them. Each norm lies strictly between zero and $N$ for $N\ge9$, so even a block gcd equal to $N$ contains a proper gcd. The product and monic remainder trees use $\widetilde O(b+\sqrt{wb})$ operations at this scale. The implementation uses the repository's uniformly fast polynomial multiplication. [Polynomial backend and exact convolution](sqrt_prime_marks.py)

Decompose each allowed row interval into maximal aligned dyadic intervals. Write $K$ for the full column span and $L$ for the number of rows. For width $w$, a row occurs in at most two tiles, and there are $O(1+K/w)$ possible positions. Cauchy–Schwarz therefore bounds the total across that level by $O(L+\sqrt{KL})$. There are $O(\log(2K))$ levels. Here $K=O(1+\sqrt H)$ and $L\le\sqrt{H/A}$, which give the displayed bit bound.

## When the principal form misses

For short $A$ and a balanced $N$ with $(-A/N)=-1$, searching every reduced primitive positive form of discriminant $-4A$ guarantees a split. For the hidden splitting prime $r$, a root $z^2=-A\bmod r$ gives the primitive form $[r,2z,(z^2+A)/r]$. Reduction preserves representation of $r$. The algorithm enumerates reduced forms directly from $A$ without computing $z$. For $[a,b,c]$ it searches the exact ellipse $a i^2+b i j+c j^2\le\lfloor\sqrt{2N}\rfloor$ with $j>0$ and signed $i$. The axis $j=0$ is covered by a gcd of the leading coefficient. Product trees work with the same general quadratic in $Y+k$.

Reduced forms have $|b|\le a\le c$ and $a i^2+b i j+c j^2\ge(ai^2+cj^2)/2$. Enumerating the forms and searching their ellipses gives $N^{1/4+o(1)}\operatorname{poly}(A)$ bit time. This is a direct norm search for short $x$. General deterministic factoring already achieves $N^{1/5+o(1)}$. For a single prescribed form and larger $A$, the faster search requires an actual hit in that form. [Harvey and Hittmeir](https://doi.org/10.1090/mcom/3708)

At $N=19\cdot23$, $x=4$, the principal ellipse contains no hit. The reduced form $3i^2-2ij+6j^2$ gives 23 at $(1,2)$. All 950 small principal-form searches agreed with direct enumeration, and the all-form search split all 468 mixed inputs in the same cohort. On the complete miss $131101\cdot163901$ at $x=1$, block evaluations fell from 10,086 to 2,681 and scalar gcds from 4,809 to zero. The maximum polynomial degree grew from 50 to 554. The 54-bit planted case at $x=8192$ returned 83,317,541 in 128 block evaluations and 59 scalar gcds. [Reproduce the comparisons](check_pooled_conic_batch.py)
