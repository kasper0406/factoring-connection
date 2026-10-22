# GPT-6.1-sol adversarial referee report on dsr_results v1 (verbatim copy of codex/scratch/dsr_2/referee_report.txt)

1. Theorem SB and the claimed consequences.

Verdict: SB(a) CORRECT; the quantitative construction in SB(b) CORRECT; its Parseval comparison CORRECT WITH FIXES; the blanket derandomizer conclusion OVERCLAIMED.

Write M=phi(N), or simply |G| in the abstract theorem. In (a), |G\W| >= delta M, so taking the complement proves the assertion, including integer rounding. A set hitting every subset of cardinality ceil(delta M) needs at least M-ceil(delta M)+1 endpoints. This is strictly greater than (1-delta)M. The adaptive argument is valid for algorithms whose only access to S is membership queries; it is not an oracle separation for the actual arithmetic predicate.

The Bernstein calculation is valid. With v=delta M and l=ln(8M), s=sqrt(6vl), v>=24l gives s<=v/2. The denominator is at most 7v/3, so the exponent is at least 18l/7, which exceeds 2l. Each of the 2(M-1)+1 bad events has probability at most 2/(8M)^2; the union has probability less than 1/(16M). The complex error is at most sqrt(2)s. The additional size condition gives 2 delta |W| <= sqrt(12 delta M l), exactly as required. There is no Bernstein or union-bound failure here.

The Parseval minimum must use the ACTUAL density alpha=|S|/M:

  beta(S) >= sqrt(|S|(M-|S|)/(M-1)).

The hypotheses only ensure alpha in [delta/2,3delta/2]. They do not ensure alpha=delta(1+o(1)): delta M / ln M is allowed to stay bounded. Even replacing alpha by delta, the displayed upper bound divided by the nominal Parseval lower bound is sqrt(48/(1-delta))*sqrt(l), as large as 8 sqrt(l), rather than 7 sqrt(l). Using the actual density, the printed estimates safely give 11 sqrt(l). These observations invalidate the stated justification of the constant 7; they are not a counterexample to the possibility of proving a sharper existence theorem. In the intended DSR regime delta=o(1), delta M/ln M tends to infinity, so 7 sqrt(l) is eventually justified. State these additional hypotheses or change the constant/comparison.

There is another quantifier error at lines 150-151: the constructed set need not have |S|>=delta M. The size-approximation theorem does not directly contradict the exact lower-density certificate written there. Fix this using a fixed-cardinality random subset and a sampling-without-replacement concentration bound, or increase the Bernoulli mean and recompute constants. Merely deleting elements from an undersized outcome does not repair it.

The assertion at lines 145-147, repeated in the opening verdict, that EVERY listed tool has a guarantee containing arbitrary density-delta accepting sets is false:

  * Nisan/INW/Nisan-Zuckerman: the usual space/width-dependent guarantees do allow arbitrary terminal predicates. With a sufficiently mixed reference walk and error below its guaranteed acceptance probability, SB(a) applies to a single generator required to work for all those predicates.
  * HPV and the width-independent BHPP result: SINGLE accepting state is an essential hypothesis. Their guarantees do not contain all terminal subsets. For HPV's PRG, summing |S| endpoint errors requires epsilon=O(1/M). BHPP's single-accept HSG cannot be converted into a many-accept HSG merely by this summation argument: HSGs do not promise additive approximation. One can instead use acceptance probability of a successful endpoint, or use their width-dependent result; neither supplies the requested short seed. Merging accepting states destroys regularity/permutation structure. See [HPV](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.ITCS.2021.7) and [BHPP](https://williamhoza.com/research/papers/BHPP/).
  * GKM fools Fourier shapes, including each endpoint character. An arbitrary endpoint indicator is a linear combination of shapes, not generally one shape. Its error bound is epsilon A'(S), after Fourier expansion. SB(a) is not a direct lower bound for GKM's stated test class. See [GKM](https://arxiv.org/abs/1506.04350).
  * KNP's main theorem concerns a FIXED finite group, equivalently constant-width permutation programs; its group-dependent constants are not a uniform theorem for G=(Z/NZ)* as N grows. De supplies more explicit group/width dependence. Neither is accurately summarized by the universal test-class sentence. See [KNP](https://iuuk.mff.cuni.cz/~koucky/papers/bwbp.pdf) and [De](https://ieeexplore.ieee.org/document/5959811/).
  * AKS/IZ amplification retains an initial random seed. Its direct parameter substitution is inadequate here; it is not itself the asserted endpoint-support impossibility theorem.

Correct conclusion: density alone cannot certify a small universal endpoint hitting set; these particular black-box applications do not prove DSR. Incorrect conclusion: generic derandomization, including every listed theorem or every possible application of it, has been refuted. IWD in the same table is itself a counterexample to the unqualified wording.

2. Proposition RW.

Verdict: CORRECT WITH FIXES concerning ranges, exact randomness, and endpoint heights.

The unweighted POSITIVE-prime measure is legitimate. One must not deduce its complex-character modulus bound from the symmetrized adjacency eigenvalue alone. JMV's proof uses the stronger estimate

  sum_{p<=y} chi(p) log p = O(sqrt(y) log y log(Ny));

partial summation gives |sum chi(p)|=O(sqrt(y) log(Ny)), including induced characters. Division by pi(y) yields the stated C log y log(Ny)/sqrt(y). Thus the measure/weight issue is repairable and does not invalidate RW. See [JMV, proof of Theorem 1.1](https://arxiv.org/pdf/0811.0647).

The deliberately loose 10^6 bound is also supportable, but is a derived bound, not JMV's printed constant. The explicit-formula remainder, truncated zero sum, and zero/pole terms in Ernvall-Hytönen-Palojärvi are safely bounded by 10^4 sqrt(x) ln^3(Nx) for psi; allowing 2*10^4 after prime powers, partial summation multiplies by at most 1/ln2+2/(ln2)^2. This remains below 10^6. Supply this derivation and its x>=2/nonprincipal-character conditions instead of calling the constant unaudited. Their final AP theorem's x>=q restriction is not the theorem to apply when y=N^{o(1)}; use their character-level explicit formula. See [explicit estimates, Theorem 29 and Lemmas 17, 21, 30](https://arxiv.org/html/2003.01925).

Fourier inversion and Cauchy-Schwarz give exactly the displayed error lambda^t sqrt(M delta(1-delta)). To make it <=delta/2 it suffices that

  t >= ln(4M(1-delta)/delta)/(2 ln(1/lambda)).

The printed t_hit omits 1-delta, making it safely larger. The pointwise formula follows from |mu^{*t}(x)-1/M| <=(M-1)lambda^t/M <=lambda^t. These lengths require lambda<1. At a small y, lambda can equal 1 and there is no finite mixing conclusion. The displayed asymptotic for t_hit also assumes delta=N^{-o(1)}, as for S_can; it is false for arbitrary densities such as delta=1/M.

For S_can and fixed b, t_hit=(1+o(1)) ln N/ln B; t_infinity=(2+o(1)) ln N/ln B. Their bit counts are the corresponding entropy/seed-space counts. Literal uniform sampling from pi(B) primes is not exactly log2 pi(B) fair bits per step. A power-of-two prime prefix, with the same asymptotic spectral bound, makes the fixed-bit claim precise. A randomized implementation must also replace unknown phi and lambda by certified bounds, for example phi<=N.

The product bound is v<=N^{1+o(1)} for t_hit, and v<=N^{2+o(1)} for t_infinity. Line 114 should distinguish these. Both are below N^3 eventually. The endpoint u<N and w=1 are valid; Omega(v)=t=o(log N). Trial division costs L^{b+o(1)} per independent trial, and O(1/delta) trials give the stated expected per-target upper bound L^{b+1/(2b)+o(1)}. None of these sufficient lengths is a lower bound on an arithmetic-specific algorithm.

3. Proposition AA and section 2.7.

Verdict: AA CORRECT; sharpness CORRECT WITH FIXES to its scope; section 2.7 OVERCLAIMED as written.

The Fourier identity and the normalization are correct. In particular

  sum_a |h(a)-delta m|^2 <= lambda^{2k} m^2 M delta(1-delta),

and division by (delta m)^2 gives |E|<=lambda^{2k}M(1-delta)/delta. Ordered tuples and repeated endpoints are allowed. Enumeration costs L^{kb+o(1)} for fixed k if only counting products; factoring/testing every endpoint by trial division adds L^{b+o(1)} per endpoint. Line 205 must not present the enumeration count as the entire constructor cost without explaining a faster batching method.

At fixed k and y=B, all UNIT residues 1<=r<N/B^k that are not B-smooth belong to E. Nonunits must be excluded. Their removal is negligible here: the number of distinct prime factors of N is O(ln N/ln B), all exceed B, and X=N/B^k is enormous. Smooth-number asymptotics give Psi(X,B)=o(X). Thus |E|>=(1-o(1))M/B^k. The ratio between the AA upper-bound exponent and this lower-bound exponent is L^{1/(2b)+o(1)}. This is a coarse comparison between two bounds, not a matching asymptotic for |E|. In particular the lower bound contains no reciprocal-density factor, and does not prove the delta^{-1} loss optimal. It is uninformative to describe AA as sharp without that qualification, especially when kb<=1/(2b) and its upper bound is vacuous.

For the q-moment, let d=pi(B) and e=|Ghat[H]|. With c_chi=M^{-1} S_hat(chi) conjugate(W_hat(chi)), the diagonal and off-diagonal estimates do yield

  sum_{q<=B}|h(q^H)-delta m|^2
    <= lambda^{2k}m^2 [d e delta(1-delta)
         + C sqrt(B) log^2(NB) A'(S)^2].

Consequently, under only A'(S)^2<=delta(1-delta)M, this PARTICULAR upper-bound method requires approximately lambda^{2k}<delta/(sqrt(B) M), up to logarithms, as well as the diagonal condition. Its N^{1+o(1)} enumeration scale is correct. Do not infer a necessary lower bound for other moments or algorithms from the size of this upper bound.

For an orbit, J must count DISTINCT residues, or the multiplicities must be retained in the estimate. Global AA implies one hit among J distinct orbit points if J>lambda^{2k}M(1-delta)/delta; replacing lambda^{2k} by about 1/m gives the quoted Jm scale. It proves one hit somewhere in the orbit, not every target, and is not a direct Parseval identity over an arbitrary short orbit.

The exceptional-set bullet confuses two goals. |E|<|T| excludes T being entirely inside E and gives at least ONE successful q. DSR requires T intersect E to be empty. Cardinality information alone ensures that only when |E|<1, irrespective of |T|. Also q -> q^H need not be injective. A B/polylog(N) lower bound on the number of distinct targets can be justified under GRH by applying Parseval to the distribution of q^H and bounding its collision probability by e/M+lambda^2; it must not be silently assumed.

4. BOX and BOX-infinity.

Verdict: BOX CORRECT; the unconditional existence assertion in BOX-infinity CORRECT WITH FIXES/clarifications; its refutation of dsr_1 section 6 CORRECT; broader DSR exclusions OVERCLAIMED if not restricted to the stated tests.

BOX(i) is simply nonwrapping multiplication retaining a prime factor >B. In (ii), (rw,1) is primitive. Any independent vector shorter than it has determinant of absolute value less than N with it, contradicting the determinant-N lattice. A dependent integral vector is an integral multiple of this primitive vector. This proves uniqueness, including ties, up to sign.

Use different symbols B_*=floor(L_{2^h}[1/2,b]) and B_N for the actual modulus-dependent threshold. With z=max(B_*,3h^2 ln^2 2)+1, sifting by primes <z excludes every prime <=max(B_N,3 ln^2 N). This remains true if B_N is unrounded, since B_N<L_{2^h}<B_*+1. N is within R of 2^h, so changing from B_* to B_N does not change the asymptotic estimates; nonsmoothness at the larger cutoff implies nonsmoothness at the smaller one. The present overloading of B obscures this valid argument.

The fundamental lemma is applicable uniformly to this TRANSLATED interval. With interval length X=R-B_*+O(1), use |A_d|=X/d+O(1), g(d)=1/d, and bounded upper/lower sieve weights supported on d<=D=R^{1/2}. Their remainder sum is O(D). The dimension-one Euler-product condition holds uniformly; s=ln D/ln z=u_h/2+o(1) tends to infinity and eventually exceeds any fixed threshold required by the lemma. Thus

  S(A,z)=X V(z)(1+O(e^{-s}))+O(D)
        =(e^{-gamma}+o(1))R/ln B_*.

No level of distribution in the location 2^h is being assumed. If the chosen textbook formulation has a divisor-weighted remainder, retain the harmless corresponding powers of log D; they are still negligible here. A transparent version of the required fundamental lemma is given in [Tao's sieve notes](https://terrytao.wordpress.com/2015/01/29/254a-supplement-5-the-linear-sieve-and-chens-theorem-optional/).

Brun-Titchmarsh really does give the needed upper bound 2R/ln R for primes in an arbitrary interval of length R, even though R=(2^h)^{o(1)}. This is an upper bound, not a short-interval prime asymptotic. Its ratio to R/ln B_* is 2/u_h=o(1). For each exponent j>=2, consecutive integer jth powers near 2^h are separated by at least a constant times 2^{h/2}; since R=2^{o(h)}, there is at most one in the interval. There are O(h) possible exponents. Removing primes and prime powers therefore leaves the main sieve term intact.

The bound Psi(R,B_*)<<R exp(-u_h/2) is valid in this range, and u_h=ceil(3 ln ln B_*) makes it O(R/(ln B_*)^{3/2}), negligible compared with R/ln B_*. A bare Euler-product Rankin estimate with a gratuitous log B factor would not suffice; quote the uniform smooth-number bound actually used, or use Hildebrand/Dickman asymptotics, valid comfortably in this range. The stated bound is recorded with a precise Tenenbaum reference in [Pollack's paper, equation (3.3)](https://www.pollack-math.net/lambdadist3.pdf).

What remains is a composite, non-prime-power, z-rough N. Since R=2^{o(h)}, it lies in (2^{h-1},2^h), H=8h exactly, r=2^h-N is nonsmooth, and r^8<N for sufficiently large h. Hence 2^H mod N=r^8. The canonical obstruction extends to integer multipliers v<=N/r^8: equality cannot occur since gcd(N,r)=1 and r>1. The stated Gauss range gives (r^8 v)^2+1<=N/2+1<N. Both conclusions follow.

For EVERY fixed k, B_N^k r^8<N eventually. Therefore the sum in dsr_1 section 6 is zero for q=2, whereas eta pi(B_N)^k delta_can is positive. This really REFUTES that missing lemma, including its quantifier k=k(b) fixed independently of N. No GRH is needed for this refutation.

It does not refute DSR. DSR allows negative integers, lifts u=r^8 wv+jN with nonzero j, |u| and |v| up to N^3, and a separate w with height as large as B_N^{log2 N}. BOX concerns a short POSITIVE combined multiplier with the canonical-positive-residue test, or the specified shortest-vector test. Signed representatives and lifted tests require separate treatment. A large w can provide wrapping even when v is small. The actual DD1 construction in DZ asks for the q=2 row, but factoring and weaker PIN designs need not literally construct that row. Likewise the orbit can wrap after O(log N) steps; BOX does not prove that such a wrap produces a hit. Replace any unconditional requirement of high multipliers for all DSR constructions by this precise test-dependent alternative.

5. Theorem FN.

Verdict: CORRECT WITH FIXES. Its coarse exponents and dense-minorant extension survive. The novelty and certificate rhetoric are OVERCLAIMED.

The displayed Fourier identity has the correct conjugations and normalization. If a translate a supp(mu^{*k}) misses S, then

  delta <= lambda^k A'(S).

There is a real gap in the prescribed choice a=ell: the least prime above B need not be a UNIT. Admissibility only rules out prime divisors <=B. A concrete admissible example at b=0.8 is

  N=9250794819020799837774052793
   =489733 * 18889465931478580854821,
  floor(L_N[1/2,0.8])=489711,
  ell=489733,
  3 ln^2 N approximately 12439.96.

Thus the printed translate is outside G. The scratch check includes a Pocklington certificate for the larger prime. Fix: choose a prime ell in (B,2B] not dividing N. Such a prime exists for every sufficiently large N because that interval has >>B/ln B primes and N has only O(ln N/ln B) prime factors. All later estimates remain valid. This same unit qualification is needed wherever the ell target is used uniformly in N.

For S_can, the printed k satisfies 2B^{k+1}<=N; equality is impossible for the relevant odd admissible N. Subtracting one from k avoids all endpoint discussion. For S_Gauss the smaller k gives ell B^k<=sqrt(N)/2, so uniqueness applies.

The logarithmic factors raised to k do NOT destroy the powers of N. With n=ln N, ln B=b sqrt(n ln n), k=O_b(sqrt(n/ln n)),

  ln(lambda^{-k}) >= k ln B/2-O(k ln n)
                   >= theta n-O_b(sqrt(n ln n)),

where theta=1/2 or 1/4. Thus delta lambda^{-k}>=N^{theta-o(1)} when delta=N^{-o(1)}. The error is generally an L-scale loss, not L^{o(1)}; it must not be discarded when claiming a precise L exponent. The note's coarse N exponents are fine.

delta_can=L^{-1/(2b)+o(1)} is correct. The assertion that coprime smooth pairs have a constant proportion of all smooth pairs requires an argument or citation; smoothness alone is not coprimality. The weaker density actually needed has a simple elementary proof. Take X=sqrt(N)/2 and t=floor(ln X/ln B). Form squarefree t-prime products using primes in (B/2,2B/3) for u and in (3B/4,B) for v. Each family has binomial(c B/ln B,t)=N^{1/2-o(1)} elements, their supports are disjoint, and all pairs are coprime and satisfy u^2+v^2<N. Their ratios are distinct by the determinant argument. Hence delta_Gauss>=N^{-o(1)}. This repairs FN without assuming the stronger displayed pair asymptotic.

For a REAL minorant f<=1_S with mean >=delta/2, expectation on the missing translate is <=0, even if f is negative outside S. Therefore

  A'(f)>= (delta/2) lambda^{-k}.

The dense-minorant extension is correct; nonnegativity of f is unnecessary. Line 290 saying minorants remain open contradicts section 2.4(i) and must be removed.

The canonical bound is sharp only at the coarse scale N^{o(1)}, since A(S)<=sqrt(|S|)<=sqrt(N). This does not prove random-set constants or norm asymptotics. The Gauss lower bound has no matching N^{1/4+o(1)} upper bound here and is not shown sharp.

For the canonical shortest-vector meaning of zroute's smooth-split set, FN does settle the speculative assertion A(S) is not L^{O(1)}, conditionally on GRH and after these repairs. It does not prove that every Bach-inspired arithmetic argument is impossible: it excludes the specified uniform small-bias/Fourier-algebra-norm certificate and its dense minorants on these sets. It says nothing comparable about every lifted/alternate test. The analytic mechanism is the elementary dual of an empty translate against a biased walk; the hard input is the already-known GRH prime bound. It should be presented as a short corollary, not as new harmonic-analytic machinery. The exact statement's priority remains unverified; relevant prior work is missing, as detailed in item 10.

6. Lemma PT and Proposition EQ-dom.

Verdict: PT's transfer CORRECT WITH FIXES to its complexity statement; its claimed implications for factoring OVERCLAIMED. EQ-dom CORRECT WITH FIXES to wording.

The proposed choice of r works for admissibility. Put d=max(b,C), A=max(L_{p^2}[alpha,d],3 ln^2(p^2)), and take r the least prime>A. For fixed 0<alpha<1, r<2A=p^{o(1)}<p eventually. Thus N=pr<p^2 has two distinct odd prime factors, each exceeding the required N-dependent cutoffs. H'=8 ceil(log2 N) is eventually in [8 log2 p,9 log2 p]; N^3=p^{3+o(1)} and the factor-base constant is b+o(1). One must define the transferred cutoff explicitly as B_N, rather than promise the exact unperturbed cutoff L_p[alpha,b].

The time L_p[alpha,c+o(1)] does not follow. Finding r has an additional L-scale cost depending on C and b. A straightforward sieve up to 2A costs L_p[alpha,2^alpha d+o(1)]; the note's naive trial-division search is no better justified at constant c. A single DSR call then costs L_p[alpha,c+o(1)]. Producing outputs for ALL pi(B_N) primes, as the statement literally says, requires that many calls, giving constant b+c, in addition to preprocessing. State a per-input-q transfer with a separate preprocessing term, or the corresponding all-q total. One can improve the preprocessing constant by choosing a threshold L_p[alpha,d+eta] and proving it dominates the eventual N cutoff, but cannot silently erase it.

The recursion observation is valid: DZ's preprocessing does not remove the chosen r. The Pollard-Strassen exponent is wrong: searching for factors <=X costs X^{1/2+o(1)} times log-size factors, so X=N^epsilon gives N^{epsilon/2+o(1)}, not N^{epsilon/4+o(1)}. The familiar N^{1/4} bound comes from searching to sqrt(N). Also do not call an attributed bound the current best without checking newer deterministic factoring results. See [Harvey-Hittmeir, discussion of Pollard-Strassen](https://link.springer.com/article/10.1007/s40993-022-00387-w).

PT is a reduction, not a proof that factoring cannot help, not a strict-strength separation, and not an oracle separation. It shows that a DSR constructor must also work on a family where N's factorization is already supplied. A hypothetical factoring-based proof could contain further arithmetic ideas. Replace 'cannot help' and 'strictly stronger' by 'no implication from deterministic factoring is established.'

For EQ-dom, (a)=>(b) can simply output phi(N)e_q after factoring N; no factorization of p-1 is needed. DD1, an affordable determinant, and ERH Miller splitting prove the reverse direction. With entries of L^{O(1)} bits, only an unspecified L constant is preserved, which is exactly the stated complexity-class equivalence. Dominant rows for every factor-base prime are sufficient, not the logically minimal data DZ needs; the small-generator PIN criterion is weaker. Complexity equivalence does not establish minimality of a particular certificate format.

7. Proposition SC and section 4.

Verdict: SC WRONG as an exact equivalence. The surrounding constructor conclusions are OVERCLAIMED.

The forward greedy factorization is valid, with a precisely maximal prefix d<=N^2: if u>N^2, d>N^2/B and c=u/d<NB. It gives an integer lift index 0<=j<B. There are B such indices for integer B and ceil(B) for real B.

For integer B, the reverse implication only proves u=dc< B N^3, not u<=N^3. For real B its immediate upper bound is ceil(B)N^3. This is a genuine false equivalence, not just a missing line. An admissible counterexample is

  N=1065023=1031*1033, B=2, a=689969,
  d=2^40=1099511627776 <= N^2,
  c=2^21=2097152,
  least_residue(a d^{-1})=1032129,
  c=1032129+N, so j=1<B.

The right side holds. But dc=2^61=2305843009213693952 exceeds N^3=1208027888215167167. Every positive 2-smooth integer <=N^3 is 2^i with i<=60, and none equals a modulo N. Indeed ord_N(2)=132870, and a=2^61 mod N. Thus the left side fails. The exact checks use integer modular arithmetic. This is a counterexample to the unconditional algebraic proposition as stated; it is not claimed as an asymptotic counterexample to GRH existence.

Exact fix: on the right add d*(least_residue(a d^{-1})+jN)<=N^3. An alternative is to state a one-way compression with a ceil(B) height slack; since B=N^{o(1)}, this gives N^{3+o(1)}, which is still not the literal DSR bound.

The multipliers are not a family of cardinality N^{2-o(1)} only: they are the inverses of smooth integers d<=N^2, whose integer labels have that order of magnitude in number, but whose residue images have cardinality at most phi(N). The test consisting of B lifts has density at most

  Psi(ceil(B)N,B)/phi(N)

as a crude upper bound; writing <=B delta_can without controlling the change in smooth density is not an exact inequality. At L scale the sufficient formulation is L^{O(1)} delta_can. Lift indicators can overlap, so no density equality follows either.

It is legitimate to say SC rewrites a search as inverse-multiplier hitting. It is not a complexity equivalence establishing 'no gain': the d family is huge, the exact height couples d and c, and a new construction might exploit that coupling. Theorem J's average fixed-parameter Johnson-class result does not rule out a worst-case arithmetic DSR algorithm; that distinction is correctly made. The Boneh dyadic-factor-base comparison is a calculation for that template, not a universal lattice lower bound.

The NP-hardness discussion needs its conditional/model qualifications. Under GRH, a suitable totalization of this polynomially verifiable DSR relation is total, and NP-hardness under the usual deterministic search reductions would imply NP=coNP under that totality assumption. Do not turn this into an unconditional exclusion based on the explicit-prime-list randomized reductions in sihard. Arbitrary explicit prime lists are a different problem. Density itself is not a general reason for non-hardness.

8. Proposition IWD.

Verdict: CORRECT WITH FIXES to implementation details; the 'would not imply circuit lower bounds' sentence is OVERCLAIMED.

With one-trial success at least delta/2, ceil(8/delta) independent trials fail with probability at most exp(-4)<0.02. Their bounded computation and success predicate have size T=L^{O(1)}. The all-sufficiently-large-lengths exponential circuit-hardness assumption in E is a standard sufficient IW assumption. An O(log T)-seed generator preserving this acceptance probability with error <1/2 has a successful seed; construction and exhaustive enumeration cost T^{O(1)}. The theorem therefore gives DSR for some constant c and, through DZ, the claimed conditional factoring theorem. See [IW](https://www.math.ias.edu/~avi/PUBLICATIONS/MYPAPERS/IW97/proc.pdf).

Uniform prime sampling from the first power-of-two prefix is sound: its endpoint cutoff is comparable to B for large B, so the same GRH estimate applies. For an actual bounded circuit, replace the noncomputable exact delta in the repetition count by a valid explicit asymptotic lower bound and pay a fixed constant slack in c. Prime-table generation, trial division, and the generator's own running time belong in T.

There is no established converse from DSR to general circuit hardness, and no such barrier is shown here. But 'it is one explicit problem' does not logically prove there can be no lower-bound implication: one explicit problem can encode a rich family, as familiar derandomization problems illustrate. Say 'no such implication or barrier is known.' The AKS analogy supports caution, not the categorical negative assertion.

9. Verdict, TH-wrap, probability, experiments, and remaining overclaims.

Verdict: 'DSR remains open here' CORRECT. The purported exact/weakest remaining obstacle and several experimental diagnoses are OVERCLAIMED.

TH-wrap is not proved equivalent to DSR, and is certainly not established as the weakest hypothesis sufficient for factoring. EQ-dom already gives an alternative sufficient certificate task, equivalent in complexity to factoring, with no prescribed smoothness/hitting format. BOX only obstructs the specified short positive multiplier tests. A successful constructor may use signed/lifted representatives, a different good set, orbit labels, or a weaker PIN design. None is forced to have the asserted form.

The proposed formulation is also incomplete. Define 'height' for a product containing inverses, reduce its numerator/denominator, and specify where the Omega budget applies. Otherwise padding m by t*t^{-1} changes the alleged height without changing its residue. A certificate for q^H m does not automatically give DSR's certificate for q^H: after clearing an inverse multiplier, u or v may exceed N^3, or w may exceed its Omega budget. Conversely DSR only bounds Omega(vw) by 4 log2 N in general, not the one-log bound printed for m. Specify how each allowed test returns a legal final (u,v,w), or directly a dominant row with the required label and affordable bit length. A map depending on N alone also needs a stated uniform-test and output-cost convention. The necessary 'wrapping height' claim is meaningful only after restricting to the canonical positive test on the BOX instances.

The exact remaining task can safely be stated as the original uniform constructive DSR problem, or as the particular arithmetic translate-hitting route one elects to pursue. The almost-all exceptional-set theorem does not establish that every route must prove the proposed TH-wrap statement. It also does not imply that all q^H must avoid a particular fixed-k E_N: BOX proves they need not, and a different constructor can use a different list/test for that target.

The <=3% assessment is subjective, as disclosed, and has no mathematical calibration. PT is an elementary closure consequence; BOX refutes one formerly proposed lemma. Neither supplies a quantitative posterior or proves DSR much harder than deterministic factoring. Do not write this percentage with the appearance of a proved upper bound.

E1's registered P1a really failed. A few large Fourier coefficients and a reduced L1 norm are consistent, but the four post-hoc examples do not prove that near-trivial characters caused the L1 ratio, or that the ratio disappears asymptotically. Mean real parts 0.33-0.86 do not mean near 1 on every factor-base prime. A lambda of 0.57 is compatible with a good expander; 'not yet an expander' is not a sound diagnosis of all these cases. Simultaneous Dirichlet approximation needs an explicit dimension/range calculation; with 36 phases and about 2.5 million characters it does not automatically force arbitrarily close alignment. GRH rules out the displayed near-1 mechanism at large L scale, but does not make smooth indicators statistically random. FN's N^{1/2-o(1)} bound permits constants below random and even L-scale suppression. The prior first-moment literature in item 10 reinforces that distinction. The actual E1 factor-base size range is 5-39; 8-36 describes only the reported post-hoc subset. Also P1c's directional pass does not confirm its registered explanation in terms of SMALL-order characters; the post-hoc characters have large order.

For E2, BOX proves m*>=Psi((N-1)/ell,B)+1 when ell is a unit; it does not prove equality or that ell is a worst target universally. I recomputed ell's first hit for ALL 20 uncensored cases and it equals the recorded m* in every one. That is additional finite-data evidence, beyond the posthoc_adv code, which only compared m* with the lower bound. The ten censored cases do not establish the exact worst target or exact m*. State the theorem and the observation separately. The seven adverse pairs really are seven pairs found in a seeded balanced-semiprime scan, not an asymptotic frequency for uniform admissible moduli or the q=2 BOX-infinity construction.

Treat censored m* as >cap, not mathematical infinity. The 24-bit lower bounds are already enough to fail P2f, so its failure survives proper censoring. For random translates, Pr[a uncovered]=(1-delta)^m is exact for each target; a Gumbel mean involving ln phi additionally assumes target-event behavior that need not hold, for example for subgroup S. D_slide's small support diagnoses this design, as the note partly acknowledges; it does not refute explicit walks as a class. Gauss/canonical DSR-target comparisons also use different b and densities, so their raw factor-of-100 difference cannot all be assigned to an inherited-vector mechanism without a matched-parameter control.

10. Literature and novelty.

Verdict: OVERCLAIMED and materially incomplete, although the main older author/year attributions are broadly sound.

There is directly relevant missing prior art on the EXACT first absolute moment that equals the canonical multiplicative Fourier-algebra norm:

  * [Hardy-Xu, Helson's conjecture for smooth numbers, 2025](https://arxiv.org/abs/2511.03430), on Steinhaus random multiplicative sums over smooth numbers, including better-than-square-root first moments and substantial quantitative suppression in relevant ranges.
  * [Hardy-Xu, Character sums over smooth numbers, July 2026](https://arxiv.org/abs/2607.00592), studying phi(q)^{-1} sum_chi |sum_{n<=x,P+(n)<=y} chi(n)|. Their upper-bound saving applies, in particular, for q>=x^{1+epsilon} in its stated y range; their theorem is not a lower bound for x=q-1. Thus it neither proves FN nor contradicts its N^{1/2-o(1)} conclusion. It does contradict the impression that smooth-character L1 behavior is an unstudied random-set question.

Read also the references there to Harper's first character moments, random multiplicative functions, and the Drappeau-Granville-Shao/Shparlinski large-smooth-character-sum literature. Add the distinction between multiplicative-character L1, additive exponential-sum norms, and sup norms of twisted sums. These are not interchangeable. The papers above provide substantive context for E1, not a theorem directly explaining these full-length finite-modulus experiments.

I did not find a prior paper stating FN's exact composite-modulus GRH lower bounds or the precise BOX-infinity theorem. That is not evidence sufficient to certify novelty. FN's empty-translate duality is elementary; BOX-infinity's sieve construction is a routine combination of classical rough-integer estimates and the congruence 2^h=r mod (2^h-r). Small-offset binary moduli themselves are old: see [Solinas, Generalized Mersenne Numbers, 1999](https://cacr.uwaterloo.ca/techreports/1999/corr99-39.pdf). That paper is background for the congruence idea, not a claimed prior proof of the rough-composite/nonsmooth-offset existence theorem. PT is an elementary padding/reduction argument; 'no prior statement found' adds little mathematical novelty.

The Ankeny/Bach citation concerns subgroup escape/generation, not a general theorem that all certificates for smooth integers require a small Fourier-algebra norm. No such general impossibility attribution was verified. Bach-Sorenson's work on random factored smooth integers is also relevant algorithmic background but does not solve prescribed-residue construction; see [Bach-Sorenson, 2020](https://arxiv.org/abs/2006.07445). Harman's unconditional smooth representatives concern a much larger power-of-modulus smoothness cutoff; they do not furnish the claimed L-scale constructor. Keep that range explicit.

The prime-finding analogy is legitimate. Cite the actual task of finding a prime in [N,2N] with a uniform running-time guarantee, not just generating an arbitrarily indexed easy subsequence of primes. The source itself is [D.H.J. Polymath, Deterministic methods to find primes](https://arxiv.org/abs/1009.3956); the journal's contributor authorship explains the Tao-Croot-Helfgott attribution. Its parity result is not a prime-finding algorithm breaking the square-root bound. Do not infer a present-day best-known claim solely from a 2012 abstract.

The sentence that deterministic subexponential relation collection in F_p* is unknown needs a DEFINITION. Arbitrary full-Q-rank relation collection is trivial: output (p-1)e_q for every factor-base prime. This is already deterministic L^{b+o(1)} output work, and EQ-dom makes the same point. Deterministic sieving of chosen index-calculus/NFS candidates is also well known, usually with heuristic yield/rank analyses. The substantive unestablished task is a uniform worst-case proved bound for short prescribed-target DSR-style witnesses, or for sufficiently informative/saturated index-calculus relations, not simply 'collecting relations.' A full-rank unsaturated lattice does not compute discrete logs. Pomerance's rigorous subexponential method is randomized; the absence of a stronger deterministic theorem in the sources checked should be reported as that limited search outcome, not an impossibility or a proof that factoring cannot assist.

Required fixes, in severity order.

FATAL: Replace Prop SC's false iff by the exact height-coupled statement, or explicitly admit the B height loss. Include the counterexample above. 'Exact same hitting problem/no gain' cannot be retained from the printed proposition.

FATAL: Withdraw the proved blanket 'generic derandomizers cannot work' verdict. SB is a universal-set support bound, not an impossibility theorem for HPV/BHPP/GKM/KNP applications or arithmetic-specific algorithms. Correct the theorem-class assertions throughout the verdict and candidate cards.

FATAL: Withdraw TH-wrap's 'exact weakest obstacle' status. Define a precise sufficient restricted route with legal final witnesses, reduced multiplier height and budgets, or return to DSR/the weaker PIN or relation task. No necessity/equivalence proof is supplied.

MAJOR: Repair FN's nonunit ell, provide a Gauss-density proof, retain the raised-logarithm loss correctly, and reconcile the contradictory dense-minorant statements. Keep the surviving GRH norm theorem.

MAJOR: Correct SB(b)'s actual-density/Parseval comparison and unsupported factor 7 in the abstract finite-group statement, and repair the exact-density consequence.

MAJOR: Correct PT's prime-selection and all-q running-time constants; change the Pollard-Strassen exponent from epsilon/4 to epsilon/2 for the cited search; remove the claimed factoring separation.

MAJOR: Separate the correct AA theorem from its coarse sharpness claim and from section 2.7's one-hit/all-target confusion, orbit distinctness, and merely sufficient template thresholds. Charge smoothness testing in constructor runtimes.

MAJOR: Add Hardy-Xu/Harper first-moment prior art; restrict Bach-certificate and deterministic-relation-collection claims to precisely defined formats. No verified novelty or circuit-lower-bound nonimplication follows from the present discussion.

MAJOR: Downgrade E1's causal/asymptotic artifact diagnosis. Separate BOX's E2 lower bound from the empirically confirmed worst target in 20 uncensored cases; handle censoring and differing smoothness parameters honestly.

MINOR: Distinguish B_* from B_N in BOX-infinity; spell out sieve uniformity, short-interval upper-bound use, and the smooth-count bound. These clarifications preserve the unconditional theorem and the refutation of dsr_1 section 6.

MINOR: State RW's lambda<1 and delta=N^{-o(1)} ranges, use exact fixed-bit sampling or qualify the entropy count, distinguish N^{1+o(1)} and N^{2+o(1)} endpoint heights, and replace subjective probability typography and universal literature claims with appropriately limited statements.

Verification artifacts: audit_checks.py and audit_checks.json in the scratch directory. Repository source files were read only. The checks cover the SC counterexample, the FN nonunit-prime example with an exact primality certificate, Fourier conjugation/normalization on a small group, and ell's first hit in all 20 uncensored E2 cases. The asymptotic sieve and GRH claims were checked analytically; toy computations are not used to prove them.
