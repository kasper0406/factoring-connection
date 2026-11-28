"""Check local rotations from each queried residue to a public occupied class."""

from itertools import permutations, product
from math import gcd, isqrt

if __package__:
    from .packed_ternary_decision import brute_sphere_count, recover_packed_symmetric_factor
else:
    from packed_ternary_decision import brute_sphere_count, recover_packed_symmetric_factor


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def reflect(vector, normal, modulus):
    scale = 2 * dot(vector, normal) * pow(dot(normal, normal), -1, modulus)
    return tuple((v - scale * w) % modulus for v, w in zip(vector, normal))


def determinant(columns):
    a, b, c = columns
    return (a[0] * (b[1] * c[2] - b[2] * c[1])
            + a[1] * (b[2] * c[0] - b[0] * c[2])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))


AXES = ((1, 0, 0), (0, 1, 0), (0, 0, 1))


def local_spin_rotation(c, y, prime, modulus):
    """Columns of a proper rotation of square spinor norm modulo prime**h."""
    target = (y, 1, 0)
    assert (dot(c, c) - dot(target, target)) % modulus == 0
    assert gcd(dot(target, target), prime) == 1
    difference = tuple(a - b for a, b in zip(c, target))
    if dot(difference, difference) % prime:
        # The first reflection maps c to the target. The second fixes it.
        normals = [difference, AXES[2]]
    else:
        # Reflection in the sum maps c to minus the target. The three
        # coordinate reflections compose to minus the identity.
        total = tuple(a + b for a, b in zip(c, target))
        normals = [total, *AXES]

    spinor = 1
    for normal in normals:
        spinor = spinor * dot(normal, normal) % prime
    if pow(spinor, (prime - 1) // 2, prime) != 1:
        # In the target's perpendicular plane, the binary norm is
        # (y²+1) + b² on the vectors (1,-y,b). Its character sum is -1.
        for b in range(prime):
            normal = (1, -y, b)
            norm = dot(normal, normal) % prime
            if pow(norm, (prime - 1) // 2, prime) == prime - 1:
                normals.extend((normal, AXES[2]))
                spinor = spinor * norm % prime
                break
    assert pow(spinor, (prime - 1) // 2, prime) == 1

    def apply(vector):
        for normal in normals:
            vector = reflect(vector, normal, modulus)
        return vector

    columns = tuple(apply(axis) for axis in AXES)
    assert apply(c) == tuple(v % modulus for v in target)
    assert determinant(columns) % modulus == 1
    for i in range(3):
        for j in range(3):
            assert dot(columns[i], columns[j]) % modulus == int(i == j)
    return columns


def check_input(n, x):
    occupied = empty = 0
    example = {}
    coordinate_residues = {}

    def decision(norm, modulus, c):
        nonlocal occupied, empty
        y = isqrt(norm - 1)
        assert y * y + 1 == norm
        ell = modulus
        while ell % 7 == 0:
            ell //= 7
        q = modulus // ell
        local_spin_rotation(c, y, 7, q)
        local_spin_rotation(c, y, ell, ell)
        count = brute_sphere_count(norm, modulus, c)
        if count:
            occupied += 1
        else:
            empty += 1
        if n == 77:
            if modulus not in coordinate_residues:
                coordinate_residues[modulus] = {
                    tuple(a * s % modulus for a, s in zip(row, signs))
                    for row in permutations((y, 1, 0))
                    for signs in product((-1, 1), repeat=3)
                }
            if not count or c not in coordinate_residues[modulus]:
                example.setdefault(ell, {}).setdefault(bool(count), (c, count))
        return count > 0

    result = recover_packed_symmetric_factor(n, x, decision)
    assert occupied and empty and result.factor > 1
    print(f"N={n}, x={x}, occupied={occupied}, empty={empty}")
    if n == 77:
        ell, counts = next((ell, counts) for ell, counts in example.items() if len(counts) == 2)
        m = result.packing_scale * ell
        print(f"  same norm {result.squared_radius} and modulus {m}: {counts}")


if __name__ == "__main__":
    for pair in ((15, 0), (35, 0), (65, 4), (77, 1), (91, 0), (183, 0)):
        check_input(*pair)
