# ldbound_4: adversarial check of the LV6/HECM transfer to design (a'''_k; mu)

Read-only repository /Volumes/git/factoring-connection. Read notes/ldbound_results.md §1, §5, §6, and in
particular the new §6.1 ("Combining with agent lv6 ... transfer to (a'''_k; mu)"). Then read
notes/lv6_results.md IN FULL (Prop LV6, Lemmas C, H, E, X, SQ, §4 field events and Lemma IR, §6 Prop BS,
§7 checklist and cost assembly).

Task: be a hostile referee of the claim that lv6's proofs (written for design (a''_k): F(g0,g1) = N,
M = N^{1/(d+1)}) transfer to (a'''_k; mu): F(g0,g1) = N w with w a random prime ~ N^{mu+1/(4k)+eta0/2},
M = N^mu, constant coefficient |a0| << M + N w / M^d << M H N^{-eta0/2}, h >= max(2nu, mu-(d-2)nu,
1-d mu+nu)+eta0, theta <= 1.
- Find EVERY place in lv6 where F(g0,g1) = N, M = N^{1/(d+1)}, the size of f's coefficients, the sizes
  of f_d / disc f / index, or the field law enters, and check the §6.1 table entry (or flag a missing row).
- Check the new content rejection (content in {1, w}) and its probability, Lemma IR with v_r = N w g1^{-d},
  N | f'(m), Lemma H's N^{O(1)} bounds, Prop BS (batch smoothness) with the new sizes, transcript invariance,
  and the order of the argument (moments under the raw law, then rejection, then Theorem K).
- Does anything in lv6 use that the leading coefficient / content / resultant is coprime to w, or that
  Res(f, g1 x - g0) = +-N (e.g. index or discriminant computations, the parity column, Lemma E)?
- Is the "Consequence" paragraph of §6.1 (Theorem A''' without the LV6 and HECM hypotheses, resting only
  on published theorems + the listed internal results) a correct statement?
Verdicts CORRECT / CORRECT WITH FIXES / WRONG per item, then an overall verdict.
