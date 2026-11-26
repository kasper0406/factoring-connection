"""Input-only factorer with a single extra predicate on the whole-disk mod-8 count."""

from argparse import ArgumentParser
from math import gcd, isqrt
from pathlib import Path
import json
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "ellipse_aggregate"))
from factor import branch_from_n
from disk import AdaptiveDisk, FixedDisks


def extract(br, a, count, y_sum):
    # One mixed positive (w,y) pair per shortest antipodal pair in one mixed lattice.
    ys = set()
    if count == 1:
        ys.add(y_sum)
    elif count == 2:
        b = isqrt(4 * a * a - br.delta)
        if b * b != 4 * a * a - br.delta or not 0 < b < a:
            raise ArithmeticError("invalid tied discriminant")
        for den in (2 * a + b, 2 * a - b):
            product, rem = divmod(a * y_sum * y_sum - br.n, den)
            gap = y_sum * y_sum - 4 * product
            if rem or product <= 0 or gap < 0:
                continue
            root = isqrt(gap)
            if root * root == gap and (y_sum - root) % 2 == 0:
                ys.update(((y_sum - root) // 2, (y_sum + root) // 2))
    for y in sorted(ys):
        v = 4 * br.n * a - br.delta * y * y
        if y <= 0 or v < 0:
            continue
        w = isqrt(v)
        if w * w != v or (w - br.delta * y) & 1:
            continue
        d = gcd((br.c * y - br.ell * w) // 2, br.n)
        if 1 < d < br.n:
            return d
    raise ArithmeticError("no verified divisor from the first-shell bit counts")


def first_mixed(br, count):
    """On [a,2a) at most three mixed antipodal pairs enter a component."""
    calls = []

    def visible(B):
        calls.append(B)
        return (count(B) - isqrt(B)) % 8

    return search(br.cap, visible, calls)


def search(cap, visible, calls):
    lo, hi = 0, 1
    while True:
        n = visible(hi)
        if n:
            break
        if hi == cap:
            raise ArithmeticError("no first shell within public bound")
        lo, hi = hi, min(2 * hi, cap)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        m = visible(mid)
        if m:
            hi, n = mid, m
        else:
            lo = mid
    return hi, n, calls


def floor_sum(n, m, a, b):
    """Sum floor((a*i+b)/m), 0<=i<n, for nonnegative arguments."""
    total = 0
    while True:
        q, a = divmod(a, m)
        total += q * n * (n - 1) // 2
        q, b = divmod(b, m)
        total += q * n
        n, b = divmod(a * n + b, m)
        if n == 0:
            return total
        a, m = m, a


def parity_first(br, oracle):
    """Detect 1 or 3 orbits by total parity, and 2 by a distinguishing y bit."""
    calls = []

    def visible(B):
        calls.append(B)
        k = isqrt(B)
        if (oracle(B) - k) & 1:
            return 1
        for i in range(isqrt(4 * br.n * B // br.delta).bit_length()):
            # Bit i of t equals floor(t/2^i) mod 2. The public y values are k*ell.
            pub = floor_sum(k, 1 << i, br.ell, br.ell)
            if (oracle(B, i) - pub) & 1:
                return 2
        return 0

    return search(br.cap, visible, calls)


def decode(br, oracle, parity=False):
    a, h, calls = parity_first(br, oracle) if parity else first_mixed(br, oracle)
    high = isqrt(4 * br.n * a // br.delta)
    k = isqrt(a)
    public_y = k * br.ell if k * k == a else 0
    total = 0
    digits = []
    for i in range(high.bit_length()):
        hits = (oracle(a, i) - oracle(a - 1, i) - ((public_y >> i) & 1)) % 8
        if parity:
            hits &= 1
        if hits > h:
            raise ArithmeticError("inconsistent first shell")
        digits.append(hits)
        total += hits << i
    if parity and h == 2:
        # First-differing bit selects one of the two distinct y coordinates.
        lowbit = total & -total
        if not lowbit:
            raise ArithmeticError("distinct positive shell pairs have distinct y")
        differing = lowbit.bit_length() - 1
        selected = 0
        for i in range(high.bit_length()):
            bits = tuple(sorted({i, differing}))
            indicator = all((public_y >> t) & 1 for t in bits)
            digit = (oracle(a, bits) - oracle(a - 1, bits) - indicator) & 1
            selected |= digit << i
        d = extract(br, a, 1, selected)
    else:
        d = extract(br, a, h, total)
    return dict(
        N=br.n,
        factors=sorted((d, br.n // d)),
        a=a,
        mixed_orbits=h,
        **({"y_xor": total} if parity else {"y_sum": total}),
        bit_counts=digits,
        radial_cutoffs=calls,
    )


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("N", type=int)
    parser.add_argument("--order", choices=("lsb", "msb", "grouped"), default="lsb")
    parser.add_argument(
        "--early-mod8",
        action="store_true",
        help="Sum out coordinate bits inside the carry calculation.",
    )
    parser.add_argument(
        "--parity-masks",
        action="store_true",
        help="Use only parity, including in the first-shell search.",
    )
    parser.add_argument(
        "--project-quotient",
        action="store_true",
        help="Fix each disk cutoff and eliminate quotient bits as well.",
    )
    args = parser.parse_args()
    if args.project_quotient and (not args.parity_masks or args.order != "lsb"):
        parser.error("--project-quotient uses --parity-masks with lsb order.")
    br, d = branch_from_n(args.N)
    if br is not None and args.parity_masks and 16 * br.cap >= args.N:
        if args.N >= 512:
            parser.error("unexpected public cap under the balanced-semiprime promise")
        d = next(p for p in range(3, isqrt(args.N) + 1, 2) if args.N % p == 0)
        result = dict(
            N=args.N, factors=sorted((d, args.N // d)), route="bounded_small_input"
        )
    elif d is not None:
        result = dict(N=args.N, factors=sorted((d, args.N // d)), route="branch_gcd")
    else:
        counter = AdaptiveDisk(
            br,
            args.order,
            early=args.early_mod8 or args.parity_masks,
            parity=args.parity_masks,
        )
        if args.project_quotient:
            counter = FixedDisks(br)
        result = decode(br, counter.count, parity=args.parity_masks)
        result.update(work=counter.metrics())
    print(json.dumps(result, indent=2))
