# Nov 26 consolidation: provenance and reproduction

Base commit 2db67a9190340630afaef7a023a9bff981bdcce9.
These notes consolidate the later ChatGPT conversation, not an original agent
worktree. Recovered five earlier ZIP archives; subsequent scratch proof and
measurement files were not present in the current environment. No missing
original raw data has been synthesized. The archive hashes below identify
the originals (the full ZIPs are retained in a separate recovery bundle
provided in the conversation; this directory contains the selected
transcribed QRP record, reconstructed notes and fresh checks, not all ZIP bytes).

| Recovered source archive | Bytes | SHA-256 |
|---|---:|---|
| qrp-experiment.zip | 6247 | 7381670f53819e062261320a90c1a7b5e674c780e8a8e9ed1fd85d8e324d1f56 |
| qrp-reciprocity-conductor.zip | 9491 | e137e002e91e00ee5d174efa4cc03eaf5e163938831f451a50718a3f4b0fdc60 |
| conductor-recovery-results.zip | 11309 | f85504dafea30b80a9b69cc9da1e307af06e94dc3bdacf4144278e70e4b4a12b |
| conductor-agreement-decoding.zip | 13886 | 54505d79d20da9a352b8b9324c918ee4eccaec7f483b216ccb59338f35d2a932 |
| jacobi-half-and-fifth.zip | 12875 | 40cb646cbfda16bc119ead55bc4a2ea59d604e459b479e55de002c4afc5e5e33 |

## Fresh checks (not the old sample)

Run from repository root with standard Python 3.8+:
    python experiments/nov26_consolidation/check_identities.py

Expected output in check_results.json. This checks the quarter-count
and reflection identities, class-number interval identities by reduced
form enumeration, slopes and grid block identities, 12,697 reduced-form
decodings, and the two explicit class-square-ambiguity examples.
The ellipse checks use known factors only to construct the diagnostic
mixed-root form and compare it with brute public counts. The public
counter here enumerates points; it is not a fast evaluator.
No claim of asymptotic complexity follows from these checks.

qrp_screen.json transcribes the original archived QRP run summary; the
original analysis plan is qrp_trial/prereg.md. The complete training
script and data are in the recovered qrp-experiment.zip; the preregistered
ML experiment was not rerun during this consolidation. Reported later
sample sizes in the notes are historical conversation reports, not newly
verified raw data.

See [shared collision note](../../notes/shared_collision_nov26.md),
[QRP note](../../notes/qrp_conductor_nov26.md), and
[roadmap](../../notes/research_roadmap_nov26.md).
