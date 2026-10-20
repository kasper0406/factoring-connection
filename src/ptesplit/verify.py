"""POST-HOC (disclosed, after codex ptesplit_2 #10): exact verification of E1 hash groups.

For every nontrivial hash group, recompute the exact power-sum vectors (int64; no overflow at these sizes, max
M*R^(M-1) < 2^63) and check that all members of a group share one vector. Equal vectors always hash equally, so
only false merges are possible; this checks for them.
"""
import json
import sys

import numpy as np

from . import census


def verify_cell(M, R, ks=(1, 2, 3)):
    arr = census.multisets(M, R)
    s_list = [M - k for k in ks]
    hs = census.fibre_hashes(arr, s_list)
    out = []
    for k in ks:
        s = M - k
        order, gs = census.group_sizes(*hs[s])
        nt = gs[gs[:, 1] >= 2]
        idx = np.concatenate([order[st:st + sz] for st, sz in nt]) if len(nt) else np.zeros(0, dtype=np.int64)
        gid = np.repeat(np.arange(len(nt)), nt[:, 1]) if len(nt) else np.zeros(0, dtype=np.int64)
        bad = 0
        sub = arr[idx]
        for i in range(1, s + 1):
            p = census.power_sum(sub, i)
            first = np.zeros(len(nt), dtype=np.int64)
            if len(nt):
                starts = np.concatenate([[0], np.cumsum(nt[:, 1])[:-1]])
                first = p[starts]
                bad += int((p != first[gid]).sum())
        out.append(dict(M=M, R=R, k=k, members_checked=int(len(idx)), mismatches=bad))
    return out


if __name__ == "__main__":
    cfg = json.load(open("experiments/configs/ptesplit.json"))
    with open("experiments/raw/ptesplit_E1_verify_posthoc.jsonl", "w") as fh:
        for M, R in cfg["E1_cells"]:
            for row in verify_cell(M, R):
                fh.write(json.dumps(row) + "\n"); print(row, flush=True)
