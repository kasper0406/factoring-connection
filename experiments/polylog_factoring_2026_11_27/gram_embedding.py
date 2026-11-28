"""Embed a binary Gram matrix primitively in Z^3 using a supplied modular root."""

from fractions import Fraction
from math import gcd, isqrt, prod
from operator import index


def _gram_schmidt(gram):
    """Exact projection coefficients and squared Gram–Schmidt lengths."""
    dimension = len(gram)
    mu = [[Fraction(0) for _ in range(dimension)] for _ in range(dimension)]
    lengths = []
    for i in range(dimension):
        for j in range(i):
            correction = sum(mu[i][k] * mu[j][k] * lengths[k] for k in range(j))
            mu[i][j] = (gram[i][j] - correction) / Fraction(lengths[j])
        length = gram[i][i] - sum(mu[i][j] ** 2 * lengths[j] for j in range(i))
        if length <= 0:
            raise ValueError("the Gram matrix must be positive definite")
        lengths.append(length)
    return mu, lengths


def orthonormal_coordinates(gram):
    """Express the first two basis vectors in an oriented orthonormal basis.

    The input must be a positive integral unimodular Gram of rank 3 or 4.
    Return two integer tuples; leave the input rows intact.
    """
    dimension = len(gram)
    if dimension not in (3, 4) or any(len(row) != dimension for row in gram):
        raise ValueError("expected a square Gram matrix of rank 3 or 4")
    gram = [list(map(index, row)) for row in gram]
    if gram != [list(column) for column in zip(*gram)]:
        raise ValueError("the Gram matrix must be symmetric")
    # Gram–Schmidt checks positive definiteness and gives the determinant.
    if prod(_gram_schmidt(gram)[1]) != 1:
        raise ValueError("the Gram matrix must be unimodular")
    # The first two input basis vectors are e1,e2. Pair them with each
    # current basis vector so its integral coordinates are ready at the end.
    pairings = [row[:2] for row in gram]
    orientation = 1
    k = 1
    while k < dimension:
        for j in range(k - 1, -1, -1):
            mu, _ = _gram_schmidt(gram)
            multiple = round(mu[k][j])
            if multiple:
                gram[k] = [x - multiple * y for x, y in zip(gram[k], gram[j])]
                for row in gram:
                    row[k] -= multiple * row[j]
                pairings[k] = [x - multiple * y for x, y in zip(pairings[k], pairings[j])]
        mu, lengths = _gram_schmidt(gram)
        if lengths[k] >= (Fraction(15, 16) - mu[k][k - 1] ** 2) * lengths[k - 1]:
            k += 1
        else:
            gram[k], gram[k - 1] = gram[k - 1], gram[k]
            for row in gram:
                row[k], row[k - 1] = row[k - 1], row[k]
            pairings[k], pairings[k - 1] = pairings[k - 1], pairings[k]
            orientation = -orientation
            k = max(k - 1, 1)

    # Successive squared GS lengths grow by at least 11/16; their product
    # is one. Thus the first integral squared norm is <=(16/11)^((rank−1)/2)
    # < 2 for rank <= 4. Size reduction makes the other basis vectors
    # orthogonal to this unit vector. Its complement is again unimodular.
    if any(gram[i][j] != int(i == j) for i in range(dimension) for j in range(dimension)):
        raise ArithmeticError("LLL failed to recover the orthonormal basis")
    pairings[-1] = [orientation * x for x in pairings[-1]]
    return tuple(pair[0] for pair in pairings), tuple(pair[1] for pair in pairings)


def embed_gram(a, b, c, root):
    """Return integer triples u,v with u.u=a, u.v=b, v.v=c and primitive u×v.

    Require a positive Gram with ac−b²=n² for positive odd n, gcd(a,n)=1,
    and root²≡−a mod n. Orientation realizes the supplied root: if t is its
    lift modulo n² and f=u×v, then −bu+av+tf is divisible by n².
    """
    a, b, c, root = map(index, (a, b, c, root))
    determinant = a * c - b * b
    if a <= 0 or determinant <= 0:
        raise ValueError("the Gram matrix must be positive definite")
    n = isqrt(determinant)
    if n * n != determinant or n % 2 == 0 or gcd(a, n) != 1:
        raise ValueError("require ac−b²=n² for odd n with gcd(a,n)=1")
    root %= n
    if (root * root + a) % n:
        raise ValueError("root does not square to −a modulo n")
    correction = -((root * root + a) // n) * pow(2 * root, -1, n) % n
    lift = root + correction * n

    # In the orthogonal sum <e1,e2> + <f>, f.f=determinant. Adjoin
    # (-b e1+a e2+lift f)/determinant. Multiplying this generator by the
    # inverse of lift gives an equivalent generator with f coefficient 1/D:
    # w=(inverse*(-b e1+a e2)+f)/D. Thus e1,e2,w are a lattice basis.
    # This is Pall's codimension-one completion, §3, (19)–(21), DOI
    # 10.4153/CJM-1949-032-9, with row (0,inverse) and ambient Gram I_3.
    inverse = pow(lift, -1, determinant)
    third = (1 + a * inverse * inverse) // determinant
    gram = [[a, b, 0], [b, c, inverse], [0, inverse, third]]
    return orthonormal_coordinates(gram)
