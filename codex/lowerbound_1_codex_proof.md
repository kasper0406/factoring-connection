# A composition theorem for smooth-relation collection

**Result.** If the $k$-th smallest fresh cofactor has log-size at least $A\sqrt{\log N\,\log(k+1)}$ at the ranks needed for extraction, the model below costs at least $L_N[1/3,(4A/3)^{2/3}]$ per target in expectation. Adaptive field selection and reuse of a common rational sieve are allowed.

Put $x=\log N$, $F_a=x^a(\log x)^{1-a}$, and $F=F_{1/3}$. Batch targets satisfy $\log N_i=(1+o(1))x$ uniformly. Write $L_N[a,c]=\exp((c+o(1))F_a)$. Work counts total operations. The result concerns relation certificates for a positive fraction of the targets in expectation.

## Model assumptions

**Public information.** Use a uniform program whose setup and candidate descriptions can be replayed from the public moduli, recorded random coins, and prior answers. Implement arithmetic in charged bit operations with independent coins and a single paid smoothness-test instruction. Constants are independent of the inputs. Tables have the same provenance, with their construction charged once. A transcript verifier checks these input wires, operation costs, multiplicative links and the extraction certificate. Only certified relation extraction earns a solved target. This rule excludes supplied CRT components, factor-dependent advice and a gcd found during setup. Replay alone says nothing about the difficulty of computing a factor from $N$.

**Tests and sharing.** A *block* specifies a target, a set of distinct candidate addresses, a smoothness bound $B=e^\mu$, and a required number $h\ge1$ of successful addresses. Open it before observing its first credited answer. A credited probe tests a fresh residual cofactor, or a bounded tuple tested jointly, of total log-size $s$. Remove known divisors first. Each revelation credits at most one address within a block, and each address only once. Order batch answers successively for H. A shared answer may guide arbitrarily many later choices without a second charge. To count an answer in several blocks, register those incidences before its revelation. If $W$ is work and $Q_i$ counts registered probes for target $i$, require $\sum_i Q_i\le S W$. We allow $S=e^{o(F)}$. In a factory the rational table can have unlimited reuse because only target-specific algebraic tests need credit.

**H [heuristic].** Given the observable history and next address, a fresh success has probability at most $p_\mu(s)=\rho(s/\mu)^{1-\epsilon_x}$, uniformly with $\epsilon_x\to0^+$, where $-\log\rho(u)\sim u\log u$. This is a stochastic oracle on labelled probes. For actual integers their smoothness is determined by their values. The arithmetic interpretation models the information obtainable before a paid test. Random cofactor estimates motivate H for blind sampling. Marginal probabilities would allow one hidden bit to make an entire field successful. Testing one address would then identify productive fields. The theorem below also works with any predictable upper probabilities satisfying its cumulative bound.

**R [model assumption].** Each successful extraction requires $h$ credits in some block. A sufficient condition is that the rows used number at least the sum of the quotas of their occupied blocks, with distinct rows assigned injectively to addresses within each block. Pigeonholing then gives R. Families using the same columns can instead pool their addresses in one block. For example, private field columns motivate a quota on the order of the number of primes up to $B$. The quantitative quota assumption is $\log h\ge(1-o(1))\mu-O(1)$. R is a separate hypothesis for a fruitful square dependency or a pinned integer relation.

## Theorem and proof

**Proved in the model.** At opening, allow each block up to $\lfloor h/4\rfloor$ exceptional addresses as automatic successes. Let $\lambda_j(k)$ be a nondecreasing lower bound on the $k$-th smallest residual log-size among the remaining distinct addresses of block $j$, uniformly over its admissible choices. Missing addresses have size infinity. For a budget of $t$ target incidences put

$$
 M_j(t)=\sum_{k=1}^{t}p_{\mu_j}(\lambda_j(k)),\qquad d_j=\lceil3h_j/4\rceil .
$$

The probability of solving a specified target within that budget is at most

$$
 \min\{1,\ t\sup_j\min[1,(eM_j(t)/d_j)^{d_j}]\}.                 \tag{1}
$$

Count a block as opened at its first registered probe. Uncredited shared work can be arbitrarily interleaved.

Indeed, conditional on opening the $a$-th probed block, let $Y$ count its ordinary successes and $P$ sum their predictable probability bounds over attempted probes, including failures. For $\theta>0$, $\exp(\theta Y-(e^\theta-1)P)$ is a nonnegative supermartingale starting at one, by the conditional Bernoulli moment bound. Distinctness and sorting give $P\le M_j(t)$ on every path up to the target budget. Stop at $Y=d_j$, localizing first at finite global time. Markov's inequality with $\theta=\log(d_j/M_j(t))$ when $M_j(t)<d_j$ gives the inner bound. Sum over at most $t$ opening positions. This proves (1) even when the menu and stopping rules depend on previous failures.

**Size–entropy consequence.** Suppose the quota inequality above holds uniformly and, for some fixed $A>0$, every ordinary block profile satisfies

$$
 \lambda_j(k)\ge A\sqrt{x\log(k+1)}.                              \tag{2}
$$

A uniform $1-o(1)$ factor in (2) is harmless. For each $c<c_A$, only ranks up to $e^{cF}$ matter, and for $h>x$ only ranks $k\ge h/8$. Then for $m$ targets with expected solved fraction bounded away from zero,

$$
 \liminf_{x\to\infty}\frac{\log(\mathbb E W/m)}{F}
 \ \ge\ c_A=(4A/3)^{2/3}.
$$

Take $t=\lfloor e^{cF}\rfloor$, $0<c<c_A$. If $h\le x$, then $\mu=O(\log x)$. Every ordinary probability is at most $\exp(-\Omega(\sqrt{x}))$, so even $t^2\sup p=o(1)$. For $h>x$, minimizing $\mu+(s/\mu)\log^+(s/\mu)$ gives $T(s)=(1+o(1))\sqrt{2s\log s}$, at $\mu\sim\sqrt{s\log s/2}$ by differentiation. Thus $\log h-\log p_\mu(s)\ge(1-o(1))\sqrt{2s\log s}$. For $h/8\le k\le t$ this exceeds $(1+\delta)\log k$ for some fixed $\delta>0$. Indeed the ratio $\sqrt{2s\log s}/z$ at $s=A\sqrt{xz}$ decreases with $z=\log k$ and tends at $z=cF$ to $\sqrt{4A/(3c^{3/2})}>1$. Consequently

$$
 M_j(t)\le h/8+h\sum_{k\ge h/8} k^{-1-\delta}=h/8+o(h).
$$

The large-quota contribution to (1) is $t e^{-\Omega(x)}=o(1)$. For every target, success therefore costs more than $t$ incidences except with probability $o(1)$. Summing expectations gives $S\mathbb E W\ge t(\mathbb E[\#\text{solved}]-o(m))$. Let $c\uparrow c_A$. The argument also tolerates $o(m)$ expected extractions outside R. Targets can be dependent. More generally, reuse $S=\exp((\sigma+o(1))F)$ subtracts at most $\sigma$ from the constant.

**Connection to SS.** For one block the opening-position factor is one. If its ordinary yield can reach a constant fraction of $h$ in expectation, some dyadic rank interval $[k,2k)$ satisfies $k p_\mu(\lambda(k))\gtrsim h/\log(2t)$. Hence $C=\log k\le\log t$ obeys $C+O(\log\log(2t))\ge(1-o(1))\sqrt{2\lambda(k)\log\lambda(k)}$. This is the SS lower fixed point with its stated scale and stability conditions. For adaptively opened private blocks, the same dyadic argument applies to $M_j(t)$. If every $M_j(t)\le d_j/(2e)$ and $\min_j d_j/\log(2t)\to\infty$, (1) vanishes. For shared columns, add cumulative counts before inversion. A union of $J$ ordinary populations satisfying (2) has profile at least $A\sqrt{x\log(1+k/J)}$, since each contributes at most $\exp(s^2/(A^2x))-1$ addresses below size $s$. The union also needs one combined exception budget. One cheap row from each of many families can otherwise fill a shared base. [Theorem SS](../../../notes/suffmap_results.md)

## Arithmetic inputs and scope

| Method | How the hypotheses enter |
|---|---|
| QS / MPQS | The conventional typical residual has log-size $(1/2+o(1))x$, giving SS exponent $1/2$. MPQS needs this estimate on the pooled supply of useful rows. |
| NFS / SNFS / tower sieving | The usual two-sided size model gives $\kappa x/D+(D/2)\log k$, hence $A=\sqrt{2\kappa}$ at the relevant ranks. Generic GNFS takes $\kappa=2$ and special-form SNFS takes $\kappa=1$. The factoring towers in B3 obey its degree-adjusted ledger. [RB2](../../../notes/resbarrier_results.md) · [B1–B3](../../../notes/heurnfs_results.md) |
| MNFS / factories | Treat each private field as a block under R. Grant all rational survivors. The usual algebraic size model $x/D+(D/2)\log k$ alone gives $A=\sqrt2$, including for a target-specific fan. Coppersmith’s factory attains a heuristic 1.6386 per target after shared preprocessing. [Coppersmith](https://research.ibm.com/publications/modifications-to-the-number-field-sieve) · [Fan factories](https://storage.googleapis.com/deepmind-media/papers/Fan_Factories_Multiple_Number_Fields_in_Factorization_Factories_and_Batch_NFS.pdf) |
| Line / lattice sieves and large primes | Charge fresh locations. Use the size of the cofactor on the selected lattice. Boundedly many primes up to $B^{O(1)}$ change the Dickman penalty by $O(\log u)$ for $u=s/\mu\to\infty$, under the usual semismooth estimate. Large-modulus selectors need a supply estimate such as LQ's selection hypothesis. [Proposition LQ](../../../notes/smoothenum_results.md) |
| Descent | In a finite descent forest each nonroot eliminates one parent prime. Under T every large prime has exactly one creator and one eliminator. Distinct roots therefore have disjoint terminal leaves, counting a direct root as a leaf. Assign one leaf per root. A root quota transfers to its block of terminal probes if these have distinct addresses. Apply (1) to the pooled **residual** profile after removing each prescribed $q$. T proves the injection only. Descent-conditioned H and a terminal size/count bound over all chosen $q$'s are additional inputs. [Assumptions HD and T](../../../notes/heurnfs_results.md) |
| $\mathbf F_2$ square dependencies / $\mathbb Z$ route | The proof uses only row demand. It covers either extraction under R. PIN can require just a few pinned coordinates, so generic full-base demand requires separate justification. [Proposition PIN](../../../notes/zroute_results.md) |

**A proved arithmetic input for unions.** Take $J$ squarefree binary product forms of common degree $D\ge3$ with nonzero resultants divisible by $N$. Then $N^2$ divides each product discriminant. Lemma C2 bounds each cumulative count by $16\exp(2s/D-2x/[D(D-1)])$ plus a box error. Discard the smallest values up to each integer error allowance. If these allowances sum to at most $h/4$, the pooled ordinary profile satisfies

$$
 \lambda(k)\ge x/(D-1)+(D/2)\log(k/(16J)).
$$

Thus $\log(16J)\le\eta\log k$ at the relevant ranks, for fixed $\eta<1$, gives $A=\sqrt{2(1-\eta)}$ asymptotically there. This handles arbitrary size-biased selection inside those boxes. Special-$q$ lattices require the transformed cofactor count. [Lemma C2 and Corollary RB2](../../../notes/resbarrier_results.md)

**A private-field certificate.** For primitive points with $b>0$ in a dyadic band $H\le b<2H$, slopes are distinct and separated by $1/(4H^2)$. They lie in $I=\{u:\lvert f(u,1)\rvert\le e^s/H^D\}$, a union of at most $D+1$ intervals. The band contributes at most $D+1+4H^2\lvert I\rvert$ points. For heights $y\in[H/2,H]$, those slopes fill area $3H^2\lvert I\rvert/8$ inside $\lvert f\rvert\le e^s$. These strips are disjoint up to boundaries. Thus a box of radius $E$ contributes $O(\mathrm{area}+D\log(2E))$, including signs and primitive points at $b=0$. For a private binary form of degree $D\ge3$ with $\log\lvert\operatorname{disc}f\rvert\ge\xi x$, the area bound in C2 now gives

$$
 \lambda(k)\ge \xi x/[2(D-1)]+(D/2)(\log k-O(1))
$$

after $O(D\log(2E))$ exceptions. This verifies $A=\sqrt{\xi}$ at growing ranks if the errors fit the private quota, even on a subset selected by rational tests. Using the typical algebraic profile for fields without this certificate needs its own justification. [Lemma FT](../../../notes/subthird_results.md) · [Area bound](../../../notes/resbarrier_results.md)

**Higher rank.** For general position the rank-$r$ count gives $s\ge\sqrt{x\log k}$ at the eligible ranks under the box-error and prefactor bounds of R3. Frame–Hölder gives related certificates where its frame exists. For a reduced noncoincident KH pair of relative module rank $\ell$ over a fixed field $k$ of degree $e$, let $n_1,n_2$ be the relative field degrees. Put $r=e\ell$, $D=e(n_1+n_2)$ and $n=\max n_i$. The larger projection dimension $b\ge\ell/2$ satisfies $b\le n$. Here $\widetilde H$ sums the log embedding maxima of the generators after subtracting fixed ideal log-norms. Ideal-relative KH gives

$$
 (2\lceil2n/b\rceil+1)\widetilde H\ge x-D\log(C_kD).
$$

If the error is at most $x/2$, the multiplier is at most $14n/\ell$, so $\widetilde H\ge rx/(28D)$. One checkable subclass uses distinct coefficient vectors $a\in\mathbb Z^r$ with residual log-size at least $\theta[\widetilde H+D\log(2+2\|a\|_\infty)]$, for fixed $\theta>0$, outside the exception allowance. At most $(2E+1)^r$ vectors have radius at most $E$. Thus $\lambda(k)\ge\theta[rx/(28D)+(D/r)\log(k+1)]$, giving (2) with $A=\theta/\sqrt7$ by AM–GM. At $\theta=1/2$, a mixed model whose other blocks meet this same floor has $c_A\approx0.39894$. Extending this subclass to sparse high-rank sampling needs a count of useful exceptional norms modulo units, subspaces and multiplicative dependencies. KH alone permits norms 1 and 8 in a high-height example. [Rank-r count](../../../notes/rank3_results.md) · [Frames](../../../notes/rankdegen_results.md) · [KH and the open counting problem](../../../notes/subthird_results.md)

## Tightness at exponent one third

**Heuristic upper checks.** With efficient enumeration and tests, the usual yield equality, concentration, fruitful dependencies and quadratic sparse linear algebra, GNFS attains $c_2=(64/9)^{1/3}$ from SS. For classical MNFS write $e,\beta,a,v$ for the $F$-scaled logarithms of box side, rational base, private base and field count. Write the degree as $\delta(x/\log x)^{1/3}$. The rational and algebraic penalties are $R=1/(3\beta\delta)$, $P=(\delta e+1/\delta)/(3a)$. The required inequalities are

$$
 2e+v-R-P\ge s=\max(\beta,a+v),\qquad
 c\ge\max(2e,2e-R+v,2s).
$$

At $y=\sqrt{13}$, choose $e=\beta=((46+13y)/108)^{1/3}$, $a=\beta(y-1)/3$, $v=\beta-a$, $\delta=\beta(4y-10)/3$. Substitution gives $R=2v, P=a, s=\beta$, hence $c=2\beta=((92+26\sqrt{13})/27)^{1/3}=1.9019\ldots$. Increasing $e$ by arbitrarily small fixed slack gives surplus relations since $2-\delta/(3a)>0$. The shared first test is priced once.

**Optimality in the classical MNFS program.** Put $b=c/2$. Eliminate $v\le s-a$ and minimize $a+(\delta e+1/\delta)/(3a)$ by AM–GM. Feasibility requires $2e\ge R+2\sqrt{(\delta e+1/\delta)/3}$. After division by $e$, the right side decreases with $e,\beta$, both at most $b$. Set them to $b$ in this necessary inequality. Writing $w=1/(3b^2\delta)<2$ and squaring gives $9b^3 w(4-8w+w^2)\ge4$. The cubic factor is maximized on the feasible interval at $w=(8-2\sqrt{13})/3$ by differentiation. This is exactly $R/\beta$ at the witness above, where each inequality is tight in the limit. Thus the MNFS constant is sharp for that program. The broader model also includes SNFS at 1.5263. [Classical cost programs](../../../notes/heurnfs_results.md)

## Three weakest points

1. Conditional cofactor randomness excludes a constructor that selects fresh smooth residuals above the Dickman rate. Prop SI and a beyond-Johnson one-witness decoder would supply selected witnesses. If their known factors are stripped, the burden shifts to the residual-size count. Theorem J covers specified Johnson certificates. [Prop SI](../../../notes/suffmap_results.md) · [Theorem J](../../../notes/sihard_results.md) · [Adaptive decoding](../../../notes/crtdecode_results.md)
2. Uniform residual-size counts must survive pooling, prefiltering and ideal normalization. The explicit KH subclass makes this requirement testable. Its extension to exceptional small norms, $\mathrm{LP}_{\rm exc}^{\rm rel}$, is open.
3. Quotas and limited reuse of credited witnesses are structural assumptions. Thin dependencies, heavily shared descent DAGs and a single witness factoring many related inputs can evade them. Repeating one modulus already defeats unrestricted amortization. Factoring-equivalent primitives remain open computational questions outside this certificate model. [Scope audit](../../../notes/completeness_audit_results.md) · [Equivalent targets](../../../notes/equivtarget_results.md)
