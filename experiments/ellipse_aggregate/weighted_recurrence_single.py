"""Factor from one positive-quadrant first moment on a fixed public disk."""

import json
from argparse import ArgumentParser
from collections import Counter
from math import gcd, isqrt

from factor import branch_from_n, positive_matches


def raw_trace(br, work):
    return sum(w for w, _ in positive_matches(br, br.cap, work)) % br.n


def decode_trace(br, raw):
    n = br.n
    k = isqrt(br.cap)
    principal = br.c * k * (k + 1) // 2
    normalized = (raw - principal) * br.ell * pow(br.c, -1, n) % n
    common = gcd(normalized, n)
    if 1 < common < n:
        return common, 0
    square = normalized**2 % n
    t = isqrt(square)
    if t > 0 and t * t == square:
        divisor = gcd(normalized - t, n)
        if 1 < divisor < n:
            return divisor, 1
    raise ArithmeticError("trace did not yield a proper divisor")


def factor_trace(n):
    br, d = branch_from_n(n)
    work = Counter()
    if d is None:
        raw = raw_trace(br, work)
        d, probes = decode_trace(br, raw)
        details = {
            "A": br.cap,
            "trace": raw,
            "square_tests": probes,
            "work": dict(work),
        }
    else:
        details = {"route": "branch_gcd", "work": dict(work)}
    if not 1 < d < n or n % d:
        raise ArithmeticError("factor failed exact verification")
    return dict(N=n, factors=sorted((d, n // d)), **details)


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("N", type=int)
    print(json.dumps(factor_trace(parser.parse_args().N), indent=2))
