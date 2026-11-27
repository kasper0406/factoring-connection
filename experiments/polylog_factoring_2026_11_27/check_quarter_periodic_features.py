"""Check discrepancy of periodic features on the scaled quarter bit."""

from math import gcd

from check_quarter_fourier import quarter


def check(p, q):
    n = p*q
    hidden = 1+p*(-2*pow(p, -1, q) % q)
    signs = [quarter(n, c*hidden) for c in range(n)]
    checked = 0
    for step in range(1, min(p, 20)):
        assert gcd(step, p) == 1
        for residue in range(step):
            values = signs[residue::step]
            for lo in range(0, len(values)-p+1):
                total = sum(values[lo:lo+p])
                assert abs(total)*q <= 10*q+8*step*p
                checked += 1
    return checked


if __name__ == "__main__":
    print(sum(check(p, q) for p, q in [(11, 13), (31, 37), (101, 139)]),
          "complete progression blocks checked")
