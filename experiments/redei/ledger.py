"""E2: bit ledger for 4-ranks over a multiplier family.

Hidden assignment sigma = (e, b, mu, {x_r}) with e=(-1/p), b=(2/p), mu=(q/p), x_r=(r/p).
Everything about q follows from public Jacobi symbols. r4 of the discriminant of Q(sqrt(-kN)) is
computed symbolically from sigma (validated against the true factorisation and PARI), and we count the
assignments consistent with all observed 4-ranks.
"""
import sys, random, json, itertools, math
sys.path.insert(0, '/Volumes/git/fc-redei/experiments/redei')
from redei import kron, fund_disc, r4_redei, rank_f2, two_sylow, rank_2k
from sympy import randprime, primerange, jacobi_symbol


def b01(s):
    return 0 if s == 1 else 1


def r4_symbolic(N, k_primes, has2, sig):
    """4-rank of Cl(fund disc of Q(sqrt(-kN))), k = prod(k_primes) * 2^has2, from sigma."""
    e, b, mu, x = sig
    m = N * (2 if has2 else 1)
    for r in k_primes:
        m *= r
    D = fund_disc(m)
    # build prime-discriminant list: labels 'p','q', public odd r, and 2-part
    labels = ['p', 'q'] + list(k_primes)
    d2 = D
    star = {}
    for r in k_primes:
        star[r] = r if r % 4 == 1 else -r
    eq = e * jacobi_symbol(-1, N)
    star_sign_p = 1 if e == 1 else -1  # p* = p if p=1 mod 4
    star_sign_q = 1 if eq == 1 else -1
    prod = star_sign_p * star_sign_q * N
    for r in k_primes:
        prod *= star[r]
    two = D // prod if D % prod == 0 else None
    assert two in (1, -4, 8, -8), (D, prod, two)
    if two != 1:
        labels.append(2)
    bq = b * jacobi_symbol(2, N)
    pq_sign = -1 if (e == -1 and eq == -1) else 1  # (p/q)(q/p)
    lam = mu * pq_sign  # (p/q)

    def chi_p(d):  # (d/p) for public d in {-1,2,r}
        if d == -1: return e
        if d == 2: return b
        return x[d]

    def chi_q(d):
        return chi_p(d) * jacobi_symbol(d, N)

    def pub_char(d, f):  # (d/.) multiplicative extension
        s = 1
        if d < 0:
            s *= f(-1); d = -d
        while d % 2 == 0:
            s *= f(2); d //= 2
        if d > 1:
            # d is a public prime here
            s *= f(d)
        return s

    def disc(lab):
        if lab == 'p': return None
        if lab == 'q': return None
        if lab == 2: return two
        return star[lab]

    def entry(li, lj):
        """Kronecker (d_j / l_i) as +-1."""
        if li in ('p', 'q'):
            chi = chi_p if li == 'p' else chi_q
            if lj == 'p':  # (p*/q)
                return (star_sign_p == -1 and eq == -1 and -1 or 1) * lam
            if lj == 'q':  # (q*/p)
                return (star_sign_q == -1 and e == -1 and -1 or 1) * mu
            return pub_char(disc(lj), chi)
        if li == 2:
            if lj == 'p': return b
            if lj == 'q': return bq
            return kron(disc(lj), 2)
        # li public odd prime r
        if lj == 'p': return x[li]
        if lj == 'q': return x[li] * jacobi_symbol(li, N)
        return kron(disc(lj), li)

    t = len(labels)
    R = [[0] * t for _ in range(t)]
    for i in range(t):
        for j in range(t):
            if i != j:
                R[i][j] = b01(entry(labels[i], labels[j]))
        R[i][i] = sum(R[i][j] for j in range(t) if j != i) % 2
    return t - 1 - rank_f2(R), D


def true_sigma(p, q, rs):
    return (int(jacobi_symbol(-1, p)), int(jacobi_symbol(2, p)), int(jacobi_symbol(q, p)),
            {r: int(jacobi_symbol(r, p)) for r in rs})


def family(rs, pairs=False):
    fam = []
    for has2 in (0, 1):
        fam.append(((), has2))
        for r in rs:
            fam.append(((r,), has2))
        if pairs:
            for r, s in itertools.combinations(rs, 2):
                fam.append(((r, s), has2))
    return fam


def ledger(N, rs, obs, fam):
    """Count assignments consistent with observed 4-ranks. Exploits locality: a member with primes K
    depends only on (e,b,mu) and x_r for r in K."""
    total = 0
    per_b0 = {}
    for e, b, mu in itertools.product((1, -1), repeat=3):
        # singleton constraints
        allowed = {}
        ok = True
        base = [f for f in fam if len(f[0]) == 0]
        xdummy = {r: 1 for r in rs}
        for f in base:
            if r4_symbolic(N, f[0], f[1], (e, b, mu, xdummy))[0] != obs[f]:
                ok = False
        if not ok:
            continue
        for r in rs:
            allowed[r] = []
            for xr in (1, -1):
                xx = dict(xdummy); xx[r] = xr
                if all(r4_symbolic(N, f[0], f[1], (e, b, mu, xx))[0] == obs[f]
                       for f in fam if f[0] == (r,)):
                    allowed[r].append(xr)
        cnt = 1
        for r in rs:
            cnt *= len(allowed[r])
        per_b0[(e, b, mu)] = {r: allowed[r] for r in rs}
        total += cnt
    return total, per_b0


def main(B=30, ntrials=200, nbits=20, seed=2):
    random.seed(seed)
    rs_all = list(primerange(3, B + 1))
    stats = []
    nvalid = 0
    for it in range(ntrials):
        p = randprime(2**(nbits - 1), 2**nbits)
        q = randprime(2**(nbits - 1), 2**nbits)
        if p == q: continue
        N = p * q
        rs = [r for r in rs_all if N % r]
        fam = family(rs)
        sig = true_sigma(p, q, rs)
        obs = {}
        for f in fam:
            a, D = r4_symbolic(N, f[0], f[1], sig)
            obs[f] = a
            if it < 20:  # validate vs real factorisation
                assert a == r4_redei(D), (p, q, f)
                if abs(D) < 10**13 and it < 5:
                    assert a == rank_2k(two_sylow(D), 2)
                nvalid += 1
        total, per = ledger(N, rs, obs, fam)
        H = 3 + len(rs)
        # bits learned beyond the 2^3 B0 guess: per consistent B0, number of undetermined x_r
        undet = {k: sum(1 for r in rs if len(v[r]) == 2) for k, v in per.items()}
        stats.append({'N': N, 'hidden_bits': H, 'consistent': total,
                      'bits_learned': H - math.log2(total), 'nB0': len(per),
                      'undet_x_true_B0': undet.get(sig[:3]), 'nr': len(rs)})
    return stats, nvalid


if __name__ == '__main__':
    out = {}
    for B in (10, 30, 100):
        st, nv = main(B=B, ntrials=150)
        bl = [s['bits_learned'] for s in st]
        out[B] = {'trials': len(st), 'validated_members': nv, 'nr_mean': sum(s['nr'] for s in st) / len(st),
                  'hidden_bits_mean': sum(s['hidden_bits'] for s in st) / len(st),
                  'bits_learned_mean': sum(bl) / len(bl), 'bits_learned_min': min(bl),
                  'nB0_hist': {k: sum(1 for s in st if s['nB0'] == k) for k in range(1, 9)},
                  'undet_x_at_true_B0_hist': {k: sum(1 for s in st if s['undet_x_true_B0'] == k)
                                              for k in range(0, 4)}}
        print(B, json.dumps(out[B]))
    json.dump(out, open('/Volumes/git/fc-redei/experiments/redei/e2_ledger.json', 'w'), indent=1)
