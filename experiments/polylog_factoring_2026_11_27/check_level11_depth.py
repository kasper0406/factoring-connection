"""Exact full and one-prime laws of the integral level-11 frames at N=19·29."""

from collections import Counter
from fractions import Fraction
from random import Random

if __package__:
    from . import level11_marking as marking
    from .check_level11_marking import crt_line, direct_frame
else:
    import level11_marking as marking
    from check_level11_marking import crt_line, direct_frame


def unit_frames(k):
    """Obtain the four or six reference automorphisms from their identity ideals."""
    inverse = marking._inverse_splitting(1, Random(1))
    result = {}
    for base in ("principal", "nonprincipal"):
        principal, witnesses = marking._class_witnesses(1, inverse, (0, 0), base)
        assert principal == (base == "principal")
        result[base] = {direct_frame(point, principal, k) for point in witnesses}
        assert len(result[base]) == (4 if principal else 6)
    return result


def project(label, modulus, units):
    """Project an entire right coset, since reducing its least member need not give its least image."""
    target, (a, b, c, d) = label
    return target, min(((a*x + b*z) % modulus, (a*y + b*w) % modulus,
                        (c*x + d*z) % modulus, (c*y + d*w) % modulus)
                       for x, y, z, w in units[target])


def collision(counts):
    return Fraction(sum(v*v for v in counts.values()), counts.total()**2)


def tv_uniform(counts, size):
    total = counts.total()
    return 1 - Fraction(sum(min(size*v, total) for v in counts.values()), size*total)


def compare_laws(groups, full, size, prime, full_spectral):
    local = [Counter(group) for group in groups]
    total = full.total()
    length = prime + 1
    supports, to_uniform, to_full, missing, collisions = [], [], [], [], []
    equal_pairs = 0
    limit = Fraction(4*prime, length**2)  # squared normalized Hecke eigenvalue bound
    for counts in local:
        assert counts.total() == length
        support = len(counts)
        coll = collision(counts)
        to_u = tv_uniform(counts, size)
        overlap = sum(min(total*v, length*full[s]) for s, v in counts.items())
        to_f = 1 - Fraction(overlap, total*length)
        outside = 1 - Fraction(sum(full[s] for s in counts), total)
        assert Fraction(1, support) <= coll <= Fraction(1, size) + (1 - Fraction(1, size))*limit
        assert max(0, 1 - Fraction(support, size)) <= to_u
        assert 4*to_u**2 <= size*coll - 1
        assert 0 <= outside <= to_f
        uniform_mass = Fraction(support, size)
        excess = max(Fraction(), 1 - outside - uniform_mass)
        assert excess**2 <= uniform_mass*(1-uniform_mass)*(size-1)*full_spectral
        supports.append(support)
        collisions.append(coll)
        equal_pairs += sum(v*(v-1) for v in counts.values())
        to_uniform.append(to_u)
        to_full.append(to_f)
        missing.append(outside)
    return ((min(supports), max(supports)), max(to_uniform), max(to_full), min(missing),
            max(collisions), equal_pairs)


def check():
    p, q = 19, 29
    n = p*q
    highest_k = 5
    native = marking.framed_line_classifier(n, highest_k, Random(n))
    top = [[native(*crt_line(p, q, i, j)) for j in range(q+1)] for i in range(p+1)]
    units = unit_frames(highest_k)
    for k in range(2, highest_k+1):
        modulus = 1 << k
        size = 5 * modulus**3 // 16
        table = [[project(label, modulus, units) for label in row] for row in top]
        full = Counter(label for row in table for label in row)
        # Double-check the lower-level projection against the native classifier
        # on each target class. Its 600 high-level classifications are reused.
        if k < highest_k:
            lower = marking.framed_line_classifier(n, k, Random(n))
            checked = set()
            for i, row in enumerate(table):
                for j, label in enumerate(row):
                    if label[0] not in checked:
                        assert lower(*crt_line(p, q, i, j)) == label
                        checked.add(label[0])
            assert checked == set(units)
        else:
            assert table == top
        full_tv = tv_uniform(full, size)
        full_coll = collision(full)
        spectral = Fraction(16*p*q, (p+1)**2 * (q+1)**2)
        assert Fraction(1, len(full)) <= full_coll <= Fraction(1, size) + (1-Fraction(1, size))*spectral
        assert max(0, 1 - Fraction(len(full), size)) <= full_tv
        assert 4*full_tv**2 <= size*full_coll - 1
        rows = compare_laws(table, full, size, q, spectral)
        columns = compare_laws(zip(*table), full, size, p, spectral)
        if modulus**2 > 4*q:
            assert rows[0] == (q+1, q+1)
        if modulus**2 > 4*p:
            assert columns[0] == (p+1, p+1)
        histogram = Counter(full.values())
        histogram[0] = size - len(full)
        equal_pairs = sum(v*(v-1) for v in full.values())
        assert rows[-1] + columns[-1] <= equal_pairs
        yield modulus, size, (len(full), full_tv, full_coll, histogram, equal_pairs), rows, columns


if __name__ == "__main__":
    for modulus, size, (support, tv, coll, histogram, equal_pairs), rows, columns in check():
        print(f"M={modulus} sector={size} full support={support}, TV to uniform={tv}, collision={coll}, "
              f"occupancies={dict(sorted(histogram.items()))}, ordered off-diagonal equal pairs={equal_pairs}, "
              f"proper determinant gcds among these={rows[-1] + columns[-1]}")
        for axis, (supports, to_u, to_f, missing, max_coll, shared) in (("rows", rows), ("columns", columns)):
            print(f"  {axis} support={supports}, max TV to uniform={to_u}, max TV to full={to_f}, "
                  f"min full mass outside support={missing}, max collision={max_coll}, shared equal pairs={shared}")
