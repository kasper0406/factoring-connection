# Commit SHA map (history rewrite 2026-10-18)

On 2026-10-18 the git history was rewritten with filter-branch to remove one 110 MB raw file, experiments/raw/trackBD_tiny.jsonl. Its gzipped copy is kept. The rewrite was needed because GitHub limits files to 100 MB.

Notes cite pre-registration commits by their ORIGINAL SHAs. This table maps old → new. The original history is preserved in ../factoring-connection-orig-history.bundle (local).

| old | new | subject |
|---|---|---|
| `de664e934571` | `5348b45f2a10` |  codex logs |
| `da54116b04f2` | `539fdcf59574` |  dcltheorem: r13 fixes — output validation (proper divisors only, N<4, band prime<N) + exhaustive N<=1e4 test; McKee-Pinch §5 as direct precursor; Carella 1308.2891; qualified sieve-barrier remark; explicit scope |
| `a204b278f13e` | `88bf85b28069` |  codex r13: fresh referee of Theorem D' (ACCEPT WITH FIXES); McKee-Pinch scans |
| `d31ac31d8790` | `bc1de326576d` |  dcltheorem: Theorem D' (deterministic almost-all N^{1/6+eps}), proof, pre-registered results, referee + critic fixes, STATUS line |
| `3550e24ebead` | `05d2bf559d3c` |  dcltheorem: pre-registration + DCL-det operational code, harness, HH baseline, tests (before pre-registered runs) |
| `700d812731be` | `0891fd07aba7` |  codex r12: DCL referee audit |
| `94deebdf8d68` | `37d1929d5b49` |  divcover: results (no construction; prime-only reformulation, separation replaces prefactoring, HS-r incidence bound), exact E1/E2 data, reviews |
| `7904079aefd0` | `a912f97f2d47` |  STATUS: DCL audit + theorem agent |
| `fabb2ef47079` | `55fb688607c4` |  multimod: results (Thm A-E, DCL scheme, pre-registered experiments, critic fixes), STATUS line |
| `81dd28ee1b13` | `a1483b23f38f` |  multimod: pre-registration + DCL code (before pre-registered runs) |
| `77d8b6708708` | `f48b36e44fec` |  divcover: pre-registration, exact GAP-cover search, structured-family drivers, tests |
| `89bf61da40c9` | `69bb2d3197f2` |  codex r11: QRH consequences for factoring |
| `dbc824e691ac` | `726003629d75` |  STATUS: phase 4 |
| `9c9c192e24ad` | `ed60bcf88174` |  user ideas: paradigms 4-6 assessment (mod p^2 formula corrected) |
| `0d72a349930c` | `9fb6efae926e` |  signsplit: STATUS line, record-count fix |
| `4424577add96` | `1468b67d424a` |  signsplit: results (REFUTED by pre-registered rule), post-hoc diagnostics, fair baseline, 2-adic shift identity |
| `091196ef502e` | `c45b32d5ce48` |  user idea: chirp quadratic exponents |
| `fc3abbb91407` | `ca06f84aaa7e` |  codex r9 tower audit (T1-T4 correct, P5 proof fixed, prior art), r10 frontier map |
| `8b981f34c1be` | `315bb5f2899b` |  user idea: Gray-code Harvey benchmark + audit |
| `a31a84f5f2d0` | `32772b17556a` |  signsplit: pre-registration, identity verification (0 failures), HH baseline + sign-split variants V1/V2/V3, tests |
| `f214f238529a` | `dc7abafd53d2` |  STATUS: phase 3 launches |
| `0028aaf8d146` | `1018147f8f4b` |  STATUS: cycloharvey done |
| `24a96e21684c` | `c1b282a87e19` |  cycloharvey: pre-registered search results, results note (verdict: refuted for this model; D1 domination + N^{1/5} lower bound) |
| `00789e720471` | `36abdbd08c20` |  codex round 8: Coppersmith-exponent bridge verification (prior art: Gao et al 2512.19076) |
| `334f17d852b3` | `a84f2ffbf525` |  signsplit: user idea + STATUS |
| `f12392ea47bb` | `c81eac922447` |  user ideas: Coppersmith-exponent bridge |
| `964fb42021a5` | `be0f6e14a2f0` |  user ideas: Schoof two-semigroup obstruction, cross-difference law |
| `b7accae35040` | `9ebeb8d949fb` |  cycloharvey: calibration gate (passed), kappa table, operational Harvey/cyclotomic searches, step-4 experiments, tests |
| `1d92a4815b92` | `d6e56568ae59` |  cycloharvey: pre-registration + calibration model, target lattices, exact covering search code |
| `cadfd10e6382` | `8b5ba6a021da` |  user idea: T=MK derivation |
| `0ddde3b5c773` | `7fbfd2efd794` |  REPORT: fold in tower128 results (s mod 128), correct order-8 claim |
| `0f76e9a06daf` | `99f56bcba5b5` |  tower128: s mod 128 via psi(p) mod 16; r8 Cl(-4N) xor r8 Cl(-8N) decides it for p=q=7 mod 8 (NEXT_STEPS item 2) |
| `d52ce6a0ac8b` | `7067c0e94902` |  pm1bsgs: evaluate p-1 stage 1 + BSGS on b^N = b^x (mod p) hybrid |
| `7c6e5ec84f41` | `567cc7ab8eb2` |  STATUS: cycloharvey agent |
| `52ec7f309103` | `40361aa36514` |  codex round 7: cyclotomic Frobenius handles for Harvey-style search; user ideas |
| `db3e79b68b46` | `ce1e4c49ff42` |  primality: analysis note, explainer, STATUS line |
| `1f0c705ce519` | `04d986cdaf76` |  epsbattery: run, audit and results of the eps/QR_2 correlation battery |
| `b02e7b43ebf5` | `69f5256b1683` |  primality: experiments (leak rates, kangaroo, iterated residues, AKS, eps certificates, battery, menu), tables, tests |
| `be7723d84a6e` | `caf656ac0fbd` |  primality: instrumented tests, exact leak formulas, kangaroo, verification; pre-register eps transcript battery |
| `779a650b7c52` | `2b891b63a279` |  epsbattery: pre-registration of the eps/QR_2 correlation battery + frozen code |
| `4064a2ad94af` | `917a3bd8c71e` |  codex round 6: eps reduction map, primality leakage literature |
| `ad7081634e8c` | `3ce8a1fe3a53` |  Phase 2 ownership |
| `1be863cca158` | `1258f95bd792` |  Final REPORT.md (revised after r5 review), status |
| `6a054028ea47` | `c2ad6dbf1d5c` |  Address r5 report review: verify items 3 and 10, rewrite NEXT_STEPS |
| `1872f22adb00` | `11c6247df928` |  Final deliverables (README, NEXT_STEPS, references, review index) and r4 verification |
| `76a3fcadd5af` | `68f1629461ee` |  codex round 4: SOS bit and class-group cards |
| `fcbcb6bf841b` | `30fdbd86b0f4` |  Status corrections after r3 audit; gitignore |
| `79d75b7c9d09` | `13391f479ae1` |  r3 audit adjudication (36 items, critic right on all) and C-2R repair-cycle results |
| `ee1d59cf43be` | `96cd69067ab3` |  Track E: partial-information decoders, information floor |
| `a8920b7604c2` | `795e68195fa9` |  C-2R repair cycle: pre-registration (arms, threshold, instances, frozen X1 extension) and search code |
| `bc9a7cb8631d` | `24aafdbdcd36` |  codex round 3 critic audit |
| `1749e1076c35` | `cc119202eb55` |  Track B/D results, codex round 2 |
| `117701087702` | `5b418c6030ff` |  Track A results: dimension-conservation theorem, section experiments |
| `a4cf8ce5a35d` | `6c3a25fb0e2f` |  Phase 1 infra, Track C results, codex round 1 notes |
