"""Rational-tangent batching inside public discriminant progressions."""

from collections import Counter
from math import isqrt

from progression_moment import classes, factor_public, fixed_cap, progression_points


def cube_root(value):
    low, high = 0, 1 << ((value.bit_length() + 2) // 3)
    while low < high:
        middle = (low + high + 1) // 2
        if middle**3 <= value:
            low = middle
        else:
            high = middle - 1
    return low


def slope(numerator, radicand, output_step, length):
    precision = 8 * length * length * (numerator + 1)
    a = numerator * precision
    b = output_step * isqrt(radicand * precision * precision)
    p0, p1, q0, q1 = 0, 1, 1, 0
    while b:
        quotient, remainder = divmod(a, b)
        p, q = quotient * p1 + p0, quotient * q1 + q0
        if q > length:
            break
        p0, p1, q0, q1 = p1, p, q1, q
        a, b = b, remainder
    return -p1, q1


def graph_points(a, c, height, x_class, z_class, work):
    """Enumerate a*x^2+c*z^2=height, with positive x,z and a*x^2<=height/2."""
    x_offset, x_step = x_class
    z_offset, z_step = z_class
    first = int(x_offset == 0)
    last = (isqrt(height // (2 * a)) - x_offset) // x_step
    length = max(1, cube_root(z_step * isqrt(c * height) // (4 * a * x_step**2)))
    if length == 1:
        for k in range(first, last + 1):
            work["direct_positions"] += 1
            x = x_offset + x_step * k
            z_squared = (height - a * x * x) // c
            z = isqrt(z_squared)
            if z > 0 and (z - z_offset) % z_step == 0 and a*x*x + c*z*z == height:
                yield x, z
        return

    for start in range(first, last + 1, length):
        work["tangent_blocks"] += 1
        end = min(last, start + length - 1)
        x = x_offset + x_step * start
        residual = height - a * x * x
        u, v = slope(a * x_step * x, c * residual, z_step, length)
        floor_height = (isqrt(v*v*residual // c) - v*z_offset) // z_step
        center = floor_height - u * start
        for intercept in range(center - 2, center + 3):
            work["line_tests"] += 1
            z_constant = z_offset * v + z_step * intercept
            z_linear = z_step * u
            x_constant = v * x_offset
            x_linear = v * x_step
            quadratic = a * x_linear**2 + c * z_linear**2
            linear = a * x_constant * x_linear + c * z_constant * z_linear
            constant = a * x_constant**2 + c * z_constant**2 - height * v**2
            discriminant = linear**2 - quadratic * constant
            if discriminant < 0:
                continue
            root = isqrt(discriminant)
            if root * root != discriminant:
                continue
            for signed_root in (root, -root) if root else (0,):
                k, remainder = divmod(signed_root - linear, quadratic)
                if remainder or not start <= k <= end:
                    continue
                z, remainder = divmod(z_constant + z_linear * k, v)
                if not remainder and z > 0 and (z - z_offset) % z_step == 0:
                    yield x_offset + x_step * k, z


def tangent_points(n, delta, cap, work):
    scale, d = (2, delta // 4) if delta % 4 == 0 else (1, delta)
    factors = factor_public(d, work)
    for level in range(1, cap + 1):
        height = 4 * n * level // scale**2
        work["levels"] += 1
        step, residues = classes(height, factors, work)
        work["progressions"] += len(residues)
        for residue in residues:
            progression = (residue, step)
            for x, y in graph_points(1, d, height, progression, (0, 1), work):
                yield scale * x, y
            for y, x in graph_points(d, 1, height, (0, 1), progression, work):
                if 2 * d * y * y < height:
                    yield scale * x, y


def hybrid_points(n, delta, cap, work):
    evaluator = tangent_points if delta**5 < n**2 else progression_points
    yield from evaluator(n, delta, cap, work)


if __name__ == "__main__":
    from argparse import ArgumentParser

    cli = ArgumentParser()
    cli.add_argument("n", type=int)
    cli.add_argument("--tangent", action="store_true")
    args = cli.parse_args()
    work = Counter()
    divisor, trace = fixed_cap(args.n, work, tangent_points if args.tangent else hybrid_points)
    print({"divisor": divisor, "trace": trace, **work})
