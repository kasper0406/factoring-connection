"""Certify distinct row vectors on the entire set of actual discovery products."""

from collections import defaultdict
import hashlib
from scoring import *


def main():
    # Digests prove inequality immediately. Verify equality by byte comparison in every shared group.
    signatures = [hashlib.sha256() for _ in CAT]

    def value(i, u):
        op, a, b = SPECS[i]
        if op == 0:
            return u[i]
        if op == 1:
            return u[a] ^ u[b]
        if op == 2:
            return u[a] & u[b]
        return u[a] | u[b]

    for n in SIZES:
        u = unaries(n, "D")
        for i, digest in enumerate(signatures):
            digest.update(value(i, u))
    groups = defaultdict(list)
    for i, h in enumerate(signatures):
        groups[h.digest()].append(i)
    shared = [g for g in groups.values() if len(g) > 1]
    for n in SIZES:
        u = unaries(n, "D")
        for group in shared:
            base = value(group[0], u)
            for i in group[1:]:
                assert np.array_equal(base, value(i, u))
    report = dict(
        rows=len(SIZES) * SPLITS["D"],
        distinct_vectors=len(groups),
        catalog=len(CAT),
        multi_member_classes=len(shared),
    )
    json.dump(report, open(ROOT / "results" / "semiprime_distinct.json", "w"), indent=2)
    print(report)
    with open(ROOT / "results" / "discovery_equivalence.tsv", "w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["test", "first_identical_test_on_discovery", "class_size"])
        for group in sorted(groups.values()):
            for test in group:
                w.writerow([test, group[0], len(group)])


if __name__ == "__main__":
    main()
