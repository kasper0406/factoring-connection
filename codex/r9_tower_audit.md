You are an adversarial auditor. Read /home/knielsen/code/factoring-connection/notes/tower_results.md and the code in src/tower/ (read-only; scratch scripts only under /tmp/claude-r9/). The agent's proofs have NOT been reviewed independently; they were only brute-force checked.

1. Audit, line by line:
   - the derivation that s mod 128 depends on ψ(p) (the even quartic character mod 16) for N ≡ 1 (mod 8), and on the order-8 character for N ≡ 5 (mod 8);
   - the general lemma "s mod 2^k depends only on (N mod 2^k, p mod 2^{k−3})";
   - Theorem T3 (s mod 128 via r8 Cl(−4N) ⊕ r8 Cl(−8N) for p ≡ q ≡ 7 mod 8), including the Rédei linearity step [−r,8r',2] + [−r,4r',2] = [2,2,−r] = ψ₂(−r). Check this against Stevenhagen arXiv:1806.06250 (use web search), including the conditions under which Rédei symbols are defined and linear;
   - T4;
   - P5 (the Chebotarev argument that single Rédei symbols carry no information about ψ(p) in the QR₂-false half);
   - the isomorphism claims Cl(−12N) ≅ Cl(−3N) and Cl(−28N) ≅ Cl(−7N) on 2-parts.
   For each: CORRECT / CORRECT WITH GAP / WRONG, with an explicit gap or counterexample. Independently reproduce key numeric claims with your own exact code where feasible: brute-force class groups for small discriminants.
2. Prior art: are T3/T4-type statements, i.e. 8-ranks of Cl(−4N), Cl(−8N) for N = pq in terms of residue symbols, already in the literature? Check Rédei–Reichardt, Stevenhagen, Milović, Koymans–Pagano, and 8-rank papers for D = −4pq / −8pq.
Output Markdown with a verdict table.
