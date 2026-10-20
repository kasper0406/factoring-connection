# lv5unif_1 — blind round: precision exponents at growing degree, and the shared-rational-side trick

You are a referee-grade analytic number theorist working read-only in the repository
/Volumes/git/factoring-connection (a research project on rigorous integer factoring). Scratch files go in
your working directory only. Use web search where useful (arXiv, journals).

## Background (read these files first)
- notes/thmA_proof.md — Theorem A: a rigorous Las Vegas L[1/2, c] factoring algorithm (internal, unrefereed),
  design (a‴_k; μ): random coprime (g0,g1) ∈ [M/2,M]², random prime w ≈ N^{μ+ν+η0/2}, random c ∈ [−H,H]^d,
  F_A = F̂_{g,w} + G·R_c with G = g1 x − g0 y, Dixon split of r² mod N into products of ≤ k primes ℓ ≈ N^{1/(2k)},
  one-level special-q descent of each ℓ in a box of radius V = √ℓ·S.
- notes/pcest_results.md §§1–2, 4.1 (Theorem FM = first moment, Lemma Kern), notes/ldbound_results.md §§1–5
  (Lemma W, Theorem L2_fam, Theorem LD* = second moment), notes/lem85_results.md (Lemma 8.5′, restricted BV),
  notes/lv6_results.md §§2–4 (linear algebra, characters: Lemmas C, H, E, X, SQ, IR), notes/hypD_results.md §§2–4
  (Theorem K, Prop PZ).
- notes/recurse_results.md, especially §3 (growing degree): §3.2 Lemma S″, §3.3 size model, §3.5 inventory.
  Recursion (Theorem RA, Lemma F) removes all shape-test cost barriers; one-level descent at L[α] needs
  d ≍ D^{1−2α}, k ≍ D^{2(1−2α)} ≲ D^α (D = log N/log log N), giving α ≥ 2/5 if the descendability density
  γ ≥ (log N)^{−O(1)} with the O(1) uniform in (d, k), and α ≥ (2+τ)/(5+2τ) if log(1/γ) ≍ d^τ log log N.
- Lee–Venkatesan, arXiv 1805.08873 (JNT 187 (2018)), §5 (Thm 5.4, Lemmas 5.13, 5.18, 5.20, 5.21) and §6, §8:
  growing degree d ≍ (log N/log log N)^{1/3}, L^{o(1)} precision, first moment only.

## Questions (answer each with proofs or precise counterexamples; label proved / heuristic / open)

1. **Precision exponent τ per ingredient.** For each of: (a) Theorem FM (incl. the LV 5.13/5.18/5.20-type
   kernel/flattening steps and the bad-moduli step via Lemma 8.5′), (b) Lemma 8.5′ re-proved in the regime
   y = B = L_N[α, a], x = N^{θ} with θ ≍ D^{2α−1} → 0, special-q primes ℓ = L_N[1−α] (e.g. L[3/5] at α = 2/5),
   (c) Lemma Kern (the constant C_d in Pr_c[Φ(c)=z] ≤ C_d|Δ|/(H²‖v∧v′‖)), (d) Theorem L2_fam and its
   (log N)^{O(C0)} loss, (e) LD* (all steps, including the "≤ 12 partners", the envelopes V₋, V_max over the
   dyadic range J ∋ ℓ, and b ∈ (V/2, V]), (f) Lemma W with η0 = η0(N) → 0, (g) Theorem K / Lemma S″ split
   counting at k ≍ D^{2/5}:
   determine *exactly where* d or k enters the loss (polynomially, like d^{O(1)} — which is (log N)^{O(1)} and
   harmless — versus exponentially, like c^d, d!, C(2k,k), (log N)^{O(d)}), and whether each exponential
   loss is intrinsic or removable by a small design change (e.g. narrowing the windows of b or ℓ). Then state
   the resulting τ and the achievable α = (2+τ)/(5+2τ).
   Important: in a PZ argument γ ≥ (E N)²/E[N²], only *relative* discrepancies between the lower-bound main term
   and the upper-bound main term count as losses. Track the size of the smoothness density ϱ at the scales used
   in the lower and in the upper bound, and the cost of shifting ϱ's argument by factors such as 2^d, d^d, or B^c
   (Hildebrand–Tenenbaum: ϱ(x/z, y)/ϱ(x, y) ≈ z^{1−α(x,y)}).
2. **Shared-rational-side multi-field trick (UNVERIFIED in recurse §3.5).** Several fields F_j share one
   rational side G = g1x − g0y (same g, independent (w_j, c_j)). Does conditioning on g decouple the algebraic
   sides enough that only L^{−o(1)} algebraic precision is needed? What conditional-on-g moments are needed,
   are they provable, and does the trick change the achievable α (in particular, does it help below 2/5, or
   only with precision)?
3. **LV6 at growing degree.** Lemma X (lv6 §3.3) uses a zero-free threshold ν₁ ≍ 1/((2d)! log Δ̄) from Stark.
   Is "character primes of (2d)!·(log N)^{O(1)} bits" a valid repair, and what does it cost? Is there a
   polynomial-bit alternative (e.g. Thorner–Zaman-type unconditional Chebotarev lower bounds), with exact
   citations?
4. **Literature.** Has anyone proved a rigorous (unconditional or under standard hypotheses only for primality)
   L[α, c] bound with α < 1/2 for factoring arbitrary N? Give exact citations if so.

Give a final summary table: ingredient | where d, k enter | τ | removable? | fix. Then your verdict on the best
rigorously achievable α in this framework and the weakest point.
