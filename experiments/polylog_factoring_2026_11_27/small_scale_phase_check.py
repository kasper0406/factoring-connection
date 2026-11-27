"""Check short-scale phase unwrapping on synthetic noisy magnitudes."""

from math import gcd, prod

from mpmath import mp


def model(n, scale, hidden, family):
    modulus = n * n
    wave = mp.cospi(mp.mpf(scale * hidden % (2 * modulus)) / modulus)
    angle = mp.mpf(scale) / modulus
    if family == 'quartic':
        return mp.sqrt((5 * mp.cospi(angle) + 4 * wave)**2 + 9 * mp.sinpi(angle)**2)
    return abs(2 * (mp.cospi(4 * angle) + wave))


def coarse_wave(n, scale, magnitude, family):
    angle = mp.mpf(scale) / (n * n)
    if family == 'quartic':
        positive = mp.sqrt(magnitude**2 - 9 * mp.sinpi(angle)**2)
        return (positive - 5 * mp.cospi(angle)) / 4
    return magnitude / 2 - mp.cospi(4 * angle)


def check_phase(n, hidden, family, noise_sign):
    with mp.workprec(2 * n.bit_length() + 80):
        calls = []

        def measure(scale):
            assert 0 < scale <= 2 * n
            calls.append(scale)
            error = noise_sign * (-1)**(scale % 3) * mp.mpf('1e-7')
            # The synthetic extension lets the check exercise all angle steps,
            # including scales where a real reduction could stop with a gcd.
            return coarse_wave(n, scale, model(n, scale, hidden, family) + error, family)

        cos_theta = measure(1)
        orientation = -1 if cos_theta < mp.mpf('-0.9') else 1
        anchor = 1

        def oriented(scale):
            return orientation**scale * measure(scale)

        if orientation == -1:
            cos_theta = -cos_theta
        while cos_theta > mp.mpf('0.9'):
            anchor *= 2
            cos_theta = oriented(anchor)

        sin_theta = mp.sqrt(1 - cos_theta**2)
        assert sin_theta > mp.mpf('0.4')
        theta = mp.acos(cos_theta)
        dilation = 1
        while 2 * anchor * dilation < n:
            dilation *= 2
            scale = anchor * dilation
            offset = -1 if gcd(scale + anchor, n) == n else 1
            cos_phase = oriented(scale)
            cos_sum = oriented(scale + offset * anchor)
            sin_phase = offset * (cos_phase * cos_theta - cos_sum) / sin_theta
            angle = mp.atan2(sin_phase, cos_phase)
            turns = mp.nint((dilation * theta - angle) / (2 * mp.pi))
            theta = (angle + 2 * mp.pi * turns) / dilation

        estimate = theta * n * n / (mp.pi * anchor)
        expected = n * n - hidden if orientation == -1 else hidden
        relative_error = abs(estimate - expected) / n
        assert relative_error < mp.mpf(1) / 64
        return len(calls), max(calls) / n, float(relative_error)


if __name__ == '__main__':
    checks = []
    synthetic_components = ((101, 109), (503, 937), (101**2, 109),
                            (153247, 254119), (1000000007, 1000000009))
    for parts in synthetic_components:
        n = prod(parts)
        modulus = n * n
        left, right = (value * value for value in parts)
        idempotent = right * pow(right, -1, left)
        for family, multiplier in (('quartic', 1), ('cubic', 4)):
            if family == 'cubic' and n < 2**14:
                continue
            phase = (multiplier * (2 * idempotent - 1)) % (2 * modulus)
            hidden = min(phase, 2 * modulus - phase)
            for numerator in (hidden, n + 1, modulus - n - 1):
                for noise_sign in (-1, 1):
                    checks.append(check_phase(n, numerator, family, noise_sign))
    print(f'{len(checks)} synthetic phase intervals verified')
    print(f'Maximum number of queried scales {max(row[0] for row in checks)}')
    print(f'Largest scale divided by N {max(row[1] for row in checks):.6g}')
    print(f'Maximum phase-coordinate error divided by N {max(row[2] for row in checks):.6g}')
