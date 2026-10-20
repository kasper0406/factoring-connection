"""Python wrapper for compsearch.c (compiled on demand into a temp dir) and a pure-Python reference."""
from __future__ import annotations
import functools, json, os, subprocess, tempfile, itertools

CS = (-2, -1, 0, 1, 2)
_SRC = os.path.join(os.path.dirname(__file__), 'compsearch.c')


@functools.lru_cache(None)
def _binary():
    d = tempfile.mkdtemp(prefix='prepread_cs_')
    out = os.path.join(d, 'compsearch')
    subprocess.check_call(['cc', '-O2', '-o', out, _SRC])
    return out


def search(k, primes, timeout=None):
    r = subprocess.run([_binary(), str(k)] + [str(p) for p in primes], capture_output=True,
                       text=True, check=True, timeout=timeout)
    return json.loads(r.stdout)


def fibre(cs, p):
    """#{x in F_p : f(x) = 0}, f = (x^2 + c_k) o ... o (x^2 + c_1), cs = (c_1..c_k). Reference."""
    sq = {}
    for x in range(p):
        sq.setdefault(x * x % p, []).append(x)
    T = {0}
    for c in reversed(cs):
        T = {x for u in T for x in sq.get((u - c) % p, ())}
        if not T:
            return 0
    return len(T)


def brute_best(k, primes):
    best, arg = -1, []
    for cs in itertools.product(CS, repeat=k):
        v = min(fibre(cs, p) for p in primes)
        if v > best:
            best, arg = v, [cs]
        elif v == best:
            arg.append(cs)
    return best, arg
