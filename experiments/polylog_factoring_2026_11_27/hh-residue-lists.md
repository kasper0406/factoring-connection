# How much does a partial answer save on balanced factoring?

Let $N=pq$ with distinct odd primes $p<q<2p$, and let $m=2^k\le N$. Suppose a supplied list of $S\ge1$ distinct odd residues modulo $m$ contains $p\bmod m$. There is a **deterministic** algorithm that finds a factor in
$$
\left(S+\min\left\{N^{1/5}(S/m)^{3/5},\ N^{1/4}S/m\right\}\right)\operatorname{poly}(\log N)
$$
bit operations. Include $Nr^{-1}\bmod m$ along with each $r$ if the list is only known to contain the residue of *some* factor. The cost to obtain the list is additional. If it misses, the algorithm can still find an incidental factor. Every returned divisor is verified. [Implementation](hh_residue_lists.py) · [Single-class method](hh-low-residue.md) · [Checks](check_hh_residue_lists.py)

For example, a list of $S=m/2$ includes every odd residue. It offers only constant-factor savings inside the one-fifth-power term, and reading a huge list can cost more than factoring directly. With six additional low bits known exactly, $m/S$ grows by 64 and the collision-search term falls by $64^{3/5}\approx12.1$. This gives a cost for weaker properties that only narrow the list.

## One baby table for the whole list

Take a public bucket parameter $K\ge1$. For each admissible residue and each of $O(K)$ geometric intervals for the smaller prime, the Harvey–Hittmeir construction gives $a,b,j$ with
$$
 aq+bp=j+m^2 i,\qquad 0\le i\le 2\lambda,\qquad
 \lambda=\left\lceil4\sqrt N/(mK)^{3/2}\right\rceil
$$
on the correct residue and interval. Since $aN+b-(aq+bp)=(p-1)(aq-b)$, the baby $\beta^{m^2 i}$ and giant $\beta^{aN+b-j}$ coincide modulo $p$ for a unit $\beta$. All $S$ classes share exactly the same baby powers. An integer equality of exponents on the target weighted sum is tested by a quadratic discriminant. Local-only coincidences are found with a product polynomial and fast multipoint evaluation. The implementation reuses the existing giant generator and converts its normalized exponents back to this common stride. [Harvey–Hittmeir, Proposition 3.3 and Algorithm 4.3](https://doi.org/10.1090/mcom/3708) · [Earlier project audit for all classes](../../notes/cycloharvey_results.md)

Choose $K$ within a constant of $\max(1,(N/(m^3 S^2))^{1/5})$. This bounds the giant count $O(SK)$ and baby count $O(1+\lambda)$ by $O(S+N^{1/5}(S/m)^{3/5})$. There is no loop over all $m/2$ odd residues. For $m/S\ge N^{1/8}$ the affine algorithm instead costs $O(S+N^{1/4}S/m)$ times a polynomial in the input bits. It covers each progression in intervals of radius proportional to $N^{1/4}$ in the divided coordinate $(p-r)/m$. [Certified affine window](coupled-rank-windows.md)

### Paying for the group element

On the remaining branch, $N^{1/5}(S/m)^{3/5}\ge N^{1/8}$. Use the deterministic order-or-factor procedure with $D=8(\lfloor N^{1/4}\rfloor+1)$. Its $\widetilde O(\sqrt D)$ cost fits. A factor ends the search. Otherwise the returned $\beta$ has global order greater than $D$. Check $\gcd(\beta^{m^2 i}-1,N)$ for $1\le i\le 2\lambda+1$. A proper gcd splits $N$. If every gcd is one, the baby powers are pairwise distinct locally at both primes, as required by the standard collision algorithm. [Implemented order setup](hh-low-residue.md#an-elementary-order-or-factor-setup-for-the-implementation)

If the first nontrivial gcd equals $N$ at $i$, the powered element has exact order $i$ at both primes. Its odd part equals the odd part of the order of $\beta$ at each prime because the stride is a power of two. Raise $\beta$ to that odd part and square successively. Different local powers of two give a proper gcd. Otherwise the two local orders are the same exact integer $r>D$, which the squarings determine. Then $p\equiv1\bmod r$ and one certified affine window suffices for all sufficiently large inputs. The implementation permits a constant number of windows including at small inputs. This handles short periods without an assumption on the input or on a random base. It specializes the common-order recovery in Gao, Feng, Hu and Pan. [Algorithm 4.6 and Proposition 4.7](https://arxiv.org/html/2512.19076#S4)

## What the fixed singleton binomials buy without answering pair queries

Suppose the lowest $k$ bits of $w=p\mathbin{\&}q$ are supplied. The fixed dyadic whole-$N$ singleton tests give exactly these bits on balanced inputs. If $t$ of them are ones, including bit zero, there are at most $2^{k-t}$ possible residues of $p$ modulo $2^k$. Each known one fixes that position on *both* factors. The list can be tightened by retaining $a$ only when $a\mathbin{\&}(N/a\bmod 2^k)=w\bmod 2^k$. This uses modular inversion of an odd residue. [Dyadic bit identity](cofactor-bit-binomials.md) · [Previous entropy lower bound for the prefix list](adaptive-binomial-lift.md#how-much-does-the-fixed-singleton-menu-leave-ambiguous)

There is a tighter bound that uses the negative singleton answers as well. In a compatible state modulo $2^j$, the two possible lifts toggle bit $j$ of *both* factors. If their new bits are equal, the choices are $00$ and $11$, and the AND bit selects exactly one. If they differ, the choices are $01$ and $10$. A supplied zero permits both, and a supplied one discards the state. Thus **on the true path, a binary branch occurs precisely at the positions where the prime bits differ**.

Cap the number of such branches per path at $B$. At each depth there are at most $2^B$ live states. To see this, start with unit mass, retain it on single-child steps, split it equally on binary steps, and discard dead ends. Every surviving state has mass at least $2^{-B}$. This proves the bound without any independence assumption on the false states. The odd residue of $p$ determines its partner $N/p\bmod 2^j$ at every depth. The code builds this tree using public arithmetic and the supplied AND only. It runs $B=0,1,\ldots$, completing each nonempty list until it finds a verified divisor. It never needs to know which positions differ. On every promised input the true residue is present by budget
$$
d=\operatorname{popcount}((p\bmod m)\mathbin{\mathrm{XOR}}(q\bmod m)).
$$
Including all smaller budgets costs a polynomial factor in $k$. For $k\le\log_2 N$ the total deterministic bit work is bounded by
$$
\left(2^d+\min\{N^{1/5}(2^d/m)^{3/5},\;N^{1/4}2^d/m\}\right)\operatorname{poly}(\log N).
$$
In fact the list never exceeds $2^{k-t}$ either. The completion may find an incidental divisor before the true state reaches the list.

For distinct primes drawn from $[U,2U]$, $d\le k/2+O(k^{2/3})$ on all but an $o(1)$ fraction of pairs at either of the following choices of $k$. For independent uniform odd residues, $d$ is binomial with parameters $k-1,1/2$ and the exceptional probability is $\exp(-\Omega(k^{1/3}))$. For primes, count exceptional integer pairs in $[U,2U]^2$. Each pair of residues modulo $2^k\le U$ occurs $O((U/2^k)^2)$ times. Dividing by $\asymp U^2/\log^2 U$ prime pairs loses at most $O(\log^2 U)$. The input-dependent $k$ has at most two possible values in each sample interval. Sorting or requiring distinctness changes only constants.

| Supplied fixed singleton information | Deterministic work on a density-one subset of the balanced samples |
| --- | --- |
| Lowest $k=\lfloor(2/9)\log_2 N\rfloor$ AND bits | $N^{2/15+o(1)}$. Here $S\le N^{1/9+o(1)}$ and the shared collision term dominates. |
| Lowest $k=\lfloor(1/4)\log_2 N\rfloor$ AND bits | $N^{1/8+o(1)}$. There are at most $N^{1/8+o(1)}$ candidate residues and the affine method needs only constantly many windows per residue. |

These comparisons count the work after receiving the singleton answers. The density-one estimates concern individual inputs. There is also a separate bound on mean work that includes the tail. For uniform odd residue pairs, $\mathbb E[2^d]=(3/2)^{k-1}$ and $\mathbb E[2^{3d/5}]=((1+2^{3/5})/2)^{k-1}$. The same counting comparison for primes loses only $O(\log^2 U)$ for these nonnegative moments. Writing $X=U^2$ and taking $k=\lfloor\theta\log_2 N\rfloor$, the mean deterministic cost on this sample is at most
$$
 X^{o(1)}\left(X^{\theta\log_2(3/2)}
     +X^{1/5+\theta\log_2((1+2^{-3/5})/2)}\right).
$$
Choose $\theta=(1/5)/(\log_2(3/2)-\log_2((1+2^{-3/5})/2))\approx0.23419$ to obtain a mean bound of $X^{0.136995+o(1)}$ (exponent rounded upward). All these bounds depend on receiving the true AND bits. A native method for obtaining them remains open. The earlier repository study treated a random set of labelled bits. Here only the unlabelled AND is supplied. [Earlier branch-and-prune experiments](../../notes/trackE_candidates.md)

More generally, if $t$ specified nonunit positions in a $k$-bit prefix are predicted to be ones on the smaller prime and at most $e$ predictions are wrong, the list has at most $2^{k-1-t}\sum_{d=0}^e\binom td$ members. This follows by enumerating which predictions fail and freely assigning all unmarked nonunit positions. The same list bound applies once this accuracy promise is justified. The cost of obtaining or validating it belongs to any proposed input-only method.

## Relation to earlier bounds

The shared baby table, modular-square stride and tangent weights are from Harvey and Hittmeir. Gao, Feng, Hu and Pan's Theorem 5.3 already gives the $m^{-3/5}$ term with an additional $N^{1/6}$ setup term when **every** prime divisor lies in one given class. Their Remark 5.5 explicitly compares that promise with finding one divisor in a class. The earlier project ledger recorded their result. Here the input is balanced and allows an arbitrary *list* for one prime. The estimate accounts for its construction from a singleton transcript, the shared search, and the low-order fallback. We make no priority claim for the underlying factoring techniques. [Gao–Feng–Hu–Pan, Theorem 5.3 and Remark 5.5](https://arxiv.org/html/2512.19076#S5) · [Published paper](https://doi.org/10.1090/mcom/4188) · [Earlier project review](../../notes/specialfam_results.md)

## Reproduction

On the 69-bit product $17179869209\cdot25769803799$, using oddness alone produced 8,885 baby steps and 12,316 giant steps. Supplying the six residues $\{9,15,25,31,45,63\}$ modulo 128 produced 2,147 babies and 2,982 giants. The list was chosen to include the smaller prime's residue. Both searches returned a divisor from a modular collision. The figures count list entries. On that same product the lowest 17 AND bits, computed from the fixture factors, left just two candidates at the first successful branch budget. The affine search found a divisor in three windows at that budget. This input has unusually few differing low bits.

The checker separately tested 471 lists containing the target, 471 other lists, 157 transcripts generated from actual small factorial ratios, and 8,847 capped trees against direct enumeration. An equal-order fixture recovered the common order 1,024 and split the number in one affine window. [Reproduce](check_hh_residue_lists.py)
