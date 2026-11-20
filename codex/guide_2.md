# guide round 2: run Gemini's full feature specification

Gemini (designer/adversary) wrote an independent spec: /Volumes/git/factoring-connection/codex/guide_spec_gemini.md
(also codex/guide_gemini.out). It has a two-tier baseline (B0 analytic rough-smooth saddle point; B1 saturated size +
local-y comparator), eight global feature families G1–G8 (coordinate S-unit arithmetic; continued fractions / root
approximation; Minkowski conjugate balance, real sign chambers, covariants; sanitized shifts and Bezout/Farey
neighbours; multiplicative candidate relations; descent/parent-seed smoothness inheritance in special-q lattices;
search-history / ultra-B cosets / other-candidate digits; raw bits / MLP), a leakage checklist (incl. L1
parallel-root-line shift leak), uncertain-stratum restriction, continuous W_s and cross-side covariance targets, and
multiple-testing control.
Task: reconcile with your round-1 pipeline (this session), implement every G-family you did not already test (at least
G2, G3, G5, G6, G8 and the continuous W_s / cross-side covariance targets), apply its leakage checks (especially L1),
pre-register before running (fresh held-out moduli), and run on the higher-bound cohort where the planted control is
detectable. Report per family: held-out gain with CIs, multiple-testing-adjusted verdict, planted-control sensitivity
for that family's feature type. Also report the descent-inheritance test (G6) prominently — "do smooth parents seed
smoother children beyond size?" is the closest thing to the user's "tame the growth" idea. ≤ 1500 words; files in cwd.
