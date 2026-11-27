"""Split a composite using integer numerators of masked unit-scale mean squares."""
from fractions import Fraction
from math import gcd
from secrets import randbelow


def try_trace_split(n, numerator, draw=randbelow):
    """One trial on n coprime to 6 with at least two distinct prime factors.

    numerator(n,u) supplies n² times the masked mean squared magnitude.
    """
    if n < 5 or gcd(n, 6) != 1:
        raise ValueError("Expected n >= 5 coprime to 6")
    masks = [1]
    while len(masks) < 2*n.bit_length():
        u = 1 + draw(n-1)
        d = gcd(n, u)
        if d > 1:
            return d
        masks += [m*u % n for m in masks]
    values = [Fraction(numerator(n, m)) for m in masks]
    half = 1
    while half < len(values):
        for start in range(0, len(values), 2*half):
            for i in range(start, start+half):
                a, b = values[i], values[i+half]
                values[i], values[i+half] = a*b, a/b
        half *= 2
    power = len(values)//2
    for target in values[1:]:
        if target <= 5**power or target >= 6**power:
            continue
        lo, hi = 5, n
        while lo < hi:
            q = (lo+hi)//2
            gain = Fraction(5*q*q+4*q+1, q*q+1)
            if gain**power > target:
                lo = q+1
            else:
                hi = q
        if 1 < lo < n and n % lo == 0:
            return lo
    return None


def count_trace_numerator(n, u):
    """Exact factor-free reference count by enumeration over residues modulo n."""
    if n < 5 or gcd(n, 6*u) != 1:
        raise ValueError("Expected n >= 5 coprime to 6 and a unit mask")
    roots = [y for y in range(n) if y*(y*y-u) % n == 0]
    total = 0
    for t in range(n):
        if t*(t-1) % n:
            continue
        choices = sum(gcd(t*(y*y-u), n) for y in roots if (1-t)*y % n == 0)
        total += choices**2
    return total


def exponent_multiple_from_cubic_count(n, root_count):
    """root_count(m) counts x(x-1)^2 = 0 modulo m², for positive m."""
    if n < 2:
        raise ValueError("Expected n >= 2")
    power = 1
    product = Fraction(1)
    while power < 3*n.bit_length()+2:
        m = n**power
        product *= Fraction(root_count(m), m)
        power *= 2
    estimate = n/product
    value = (2*estimate.numerator + estimate.denominator)//(2*estimate.denominator)
    return n*value


def exponent_multiple_from_fixed_trace(n, trace):
    """trace(m) returns a rational approximation to the mask-one mean square."""
    from mpmath import mp

    if n < 5 or gcd(n, 6) != 1:
        raise ValueError("Expected n >= 5 coprime to 6")
    bits = n.bit_length()
    limit = 16*(bits+1)
    with mp.workprec(4*bits + 16*limit.bit_length() + 32):
        def read(j):
            value = Fraction(trace(n**j))
            return mp.mpf(value.numerator)/value.denominator

        last = read(limit)
        baseline = int(mp.nint(last))
        # log(1+4x/5+x²/5) has coefficient -(alpha^i+beta^i)/i.
        moments = [mp.mpf(2), -mp.mpf(4)/5]
        for j in range(2, limit+1):
            moments.append(-(4*moments[-1]+moments[-2])/5)
        coeff = [0] + [-moments[j]/j for j in range(1, limit+1)]
        residual = [0] + [-mp.mpf(1)/j for j in range(1, limit+1)]
        result = mp.mpf(0)
        for j in range(1, limit+1):
            weight = residual[j]/coeff[1]
            value = last if j == limit else read(j)
            result += weight*mp.log(value/baseline)
            for multiple in range(2*j, limit+1, j):
                residual[multiple] -= weight*coeff[multiple//j]
        return n*int(mp.nint(n*mp.exp(result)))


def fixed_trace_interval(first, last):
    """Yield (n, exact mask-one trace) for n coprime to 6 in a closed interval."""
    from math import isqrt

    if first < 1 or last < first:
        raise ValueError("Expected a nonempty interval of positive integers")
    remainder = list(range(first, last+1))
    numerator = [1]*len(remainder)
    limit = isqrt(last)
    sieve = bytearray([1])*(limit+1)
    for p in range(2, limit+1):
        if not sieve[p]:
            continue
        if p*p <= limit:
            sieve[p*p::p] = bytes((limit-p*p)//p+1)
        if p < 5:
            continue
        for i in range((-first) % p, len(remainder), p):
            if gcd(first+i, 6) != 1:
                continue
            q = p
            remainder[i] //= p
            while remainder[i] % p == 0:
                q *= p
                remainder[i] //= p
            numerator[i] *= 5*q*q + 4*q + 1
    for i, q in enumerate(remainder):
        n = first+i
        if gcd(n, 6) == 1:
            if q > 1:
                numerator[i] *= 5*q*q + 4*q + 1
            yield n, Fraction(numerator[i], n*n)


def distinct_count_from_coarse_trace(n, trace):
    """Use one mask-one trace estimate with at most 25% relative error."""
    if n < 1 or gcd(n, 6) != 1:
        raise ValueError("Expected positive n coprime to 6")
    estimate = Fraction(trace(n*n))
    count, threshold = 0, 2
    while estimate > threshold:
        count += 1
        threshold *= 5
    return count


def coarse_trace_from_distinct_count(n, count, inverse_error=4):
    """Relative error ≤1/inverse_error given the distinct-prime count of n."""
    from sympy import primerange

    if n < 1 or gcd(n, 6) != 1 or inverse_error < 1:
        raise ValueError("Expected positive n coprime to 6 and positive inverse error")
    estimate = Fraction(5**count)
    for p in primerange(5, n.bit_length()*inverse_error+1):
        if n % p:
            continue
        q = 1
        while n % p == 0:
            q *= p
            n //= p
        estimate *= Fraction(5*q*q+4*q+1, 5*q*q)
        if n == 1:
            break
    return estimate


def large_factor_hint(n, mixed_trace):
    """Invert the mixed trace assuming the smaller prime has positive mask sign."""
    from math import isqrt

    if n < 144:
        raise ValueError("Use direct division for small semiprimes")
    target = Fraction(mixed_trace)
    low, high = isqrt(n-1)+1, n
    while low < high:
        z = (low+high)//2
        numerator = (5*n*n+4*z*n+z*z)*(z*z+1)
        denominator = n*n*z*z
        if numerator*target.denominator < denominator*target.numerator:
            low = z+1
        else:
            high = z
    return low


def split_from_mixed_trace(n, estimate, absolute_error, metrics=None):
    """Try the smaller-positive orientation with a certified absolute trace error."""
    from math import isqrt
    from src.trackE.coppersmith import decode_interval

    error = Fraction(absolute_error)
    if error < 0:
        raise ValueError("Expected a nonnegative error bound")
    hint = large_factor_hint(n, estimate)
    radius = 1 + n*error//3
    lower = max(isqrt(n-1)+1, hint-radius)
    upper = min(n-1, hint+radius)
    return decode_interval(n, lower, upper-lower, metrics=metrics, enforce_rigorous=True)


def semiprime_trace_from_samples(samples, error, precision):
    """Estimate the unit-scale trace from nine consecutive absolute mean-square samples.

    samples[j-1] is |J_n(j,u)|² within absolute error; use primes greater than 9.
    precision is the working binary precision, including numerical guard bits.
    """
    from mpmath import mp

    if len(samples) != 9 or error <= 0:
        raise ValueError("Expected nine samples and a positive absolute error bound")
    with mp.workprec(precision):
        values = [mp.mpf(0), *samples]
        off = mp.matrix([[values[abs(i-j)] for j in range(10)] for i in range(10)])
        diagonal = -mp.eigsy(off, eigvals_only=True)[0]
        block = off[:9, :9] + (diagonal + 20*error)*mp.eye(9)
        ones = mp.matrix([1]*9)
        return 1/sum(mp.lu_solve(block, ones))


def split_from_five_complex_samples(n, samples):
    """Five calibrated mask-one samples at scales 1..5 on a two-component base."""
    from mpmath import mp

    if n < 7 or gcd(n, 30) != 1 or len(samples) != 5:
        raise ValueError("Expected a base coprime to 30 and five complex samples")
    modulus = n*n
    with mp.workprec(12*n.bit_length()+80):
        z = mp.expjpi(mp.mpf(2)/modulus)
        residual = [value-z**j for j, value in enumerate(samples, 1)]
        diff = [b-a for a, b in zip(residual, residual[1:])]
        middle = max((1, 2), key=lambda i: abs(diff[i]))
        if not diff[middle]:
            return None
        root_sum = (diff[middle+1]+z*diff[middle-1])/diff[middle]
        root = (root_sum + mp.sqrt(root_sum**2-4*z))/2
        exponent = int(mp.nint(mp.arg(root)*modulus/(2*mp.pi)))
        factor = gcd(exponent, n)
        return factor if 1 < factor < n else None


def split_from_one_complex_sample(n, sample):
    """Decode calibrated J_n(1,1) for two coprime prime-power components >=5."""
    from math import isqrt
    from mpmath import mp
    from sympy import Poly, Rational, symbols

    if n < 5 or gcd(n, 6) != 1:
        raise ValueError("Expected n >= 5 coprime to 6")
    bits = 4*n.bit_length()+64
    with mp.workprec(bits+64):
        h = mp.expjpi(-mp.mpf(1)/(n*n))
        rotated = h*sample-h.conjugate()

        def dyadic(value):
            return Rational(int(mp.nint(value*(1 << bits))), 1 << bits)

        u, v = dyadic(rotated.real), dyadic(rotated.imag)
        c, d = dyadic(h.real), dyadic(h.imag)
    k = symbols("k")
    s = 4+k/n
    product = 2*s-4+Rational(1, n)
    delta = (s-4)**2-Rational(4, n)
    quartic = Poly(delta*(u-product*c)**2+s**2*(v-product*d)**2-s**2*delta, k)
    tolerance = Rational(1, 10**6*n**4)
    for level in (-tolerance, tolerance):
        for (lower, _), _mult in (quartic-level).intervals(eps=Rational(1, 8)):
            center = int(lower//1)
            for total in range(max(1, center-1), center+3):
                discriminant = total*total-4*n
                if discriminant >= 0:
                    root = isqrt(discriminant)
                    factor = (total-root)//2
                    if root*root == discriminant and 1 < factor < n and n % factor == 0:
                        return factor
    return None


def fixed_sample_root_sum(n, precision):
    """Compute J_n(1,1) without supplied factors, scanning n residue classes."""
    from mpmath import mp

    with mp.workprec(precision):
        return mp.fsum(mp.expjpi(2*mp.mpf(e)/(n*n))/h for e, h in fixed_sample_terms(n))


def fixed_sample_terms(n):
    """Yield the projector and weight denominator of each admissible residue."""
    if n < 5 or gcd(n, 6) != 1:
        raise ValueError("Expected n >= 5 coprime to 6")
    for w in range(n):
        if w*(w-1)*(w+1)*(w-2) % n:
            continue
        g = gcd(w, n)
        outside = (n//g)**2
        yield outside*pow(outside, -1, g*g), gcd(w-2, n)


def fixed_lift_sample_table(n):
    """Numerical mask-one reference table indexed by k at scales 1+n*k."""
    import numpy as np

    if n < 5 or gcd(n, 6) != 1:
        raise ValueError("Expected n >= 5 coprime to 6")
    bins = np.zeros(n, dtype=complex)
    for e, h in fixed_sample_terms(n):
        bins[e % n] += np.exp(2j*np.pi*(e/(n*n)))/h
    return np.fft.ifft(bins)*n



def split_from_one_kloosterman(n, normalized_sum):
    """Decode K_(n²)(1,1)/n on odd bases with two distinct prime factors."""
    from mpmath import mp

    if n < 3 or n % 2 == 0:
        raise ValueError("Expected an odd base")
    modulus = n*n
    with mp.workprec(5*n.bit_length()+64):
        cosine = normalized_sum/2-mp.cospi(mp.mpf(4)/modulus)
        cosine = max(-1, min(1, cosine))
        index = int(mp.nint(modulus*mp.acos(cosine)/(2*mp.pi)))
        root = index*pow(2, -1, modulus) % modulus
        divisor = gcd(root-1, n)
        return divisor if 1 < divisor < n else None


def split_from_kloosterman_quotient(n, normalized_sum, error_units=1):
    """Decode with absolute error ≤error_units/n²; enumerate O(error_units) integers."""
    from math import isqrt
    from mpmath import mp

    if n < 3 or n % 2 == 0 or not isinstance(error_units, int) or error_units < 1:
        raise ValueError("Expected an odd base and a positive integer error budget")
    modulus = n*n
    with mp.workprec(4*n.bit_length()+64):
        cosine = mp.mpf(normalized_sum)/2-mp.cospi(mp.mpf(4)/modulus)
        angle = mp.acos(min(1, abs(cosine)))
        estimate = modulus*(angle/mp.pi)**2-mp.mpf(16)/modulus
        center = int(mp.nint(estimate))
    radius = (error_units+5)//6+1
    for quotient in range(max(0, center-radius), center+radius+1):
        square = 16+modulus*quotient
        root = isqrt(square)
        if root*root == square:
            divisor = gcd(root-4, n)
            if 1 < divisor < n:
                return divisor
    return None


def split_from_scaled_kloosterman(n, scale, raw_sum):
    """One trial with raw absolute error <= 1/16. Requires python-flint."""
    from mpmath import mp

    if n < 17 or n % 2 == 0 or gcd(scale, n) != 1:
        raise ValueError("Expected odd n >= 17 and a unit scale")
    modulus = n*n
    with mp.workprec(4*n.bit_length()+64):
        cosine = mp.mpf(raw_sum)/(2*n)-mp.cospi(mp.mpf(4)*(scale % modulus)/modulus)
        if abs(cosine) > mp.mpf(3)/4:
            return None
        center = int(mp.nint(modulus*mp.acos(cosine)/(2*mp.pi)))
    return _split_quadratic_interval(n, scale, center)


def _split_quadratic_interval(n, scale, center):
    from flint import fmpz_mat, fmpz_poly

    modulus = n*n
    radius = (n+63)//64+1
    m = 1
    while True:
        dim = 2*m+2
        det = modulus**(m*(m+1))*radius**(dim*(dim-1)//2)
        if 2**(dim*(dim-1))*det**4*dim**(2*dim) < modulus**(4*m*dim):
            break
        m += 1
    f = fmpz_poly([(center*center-4*scale*scale) % modulus, (2*center) % modulus, 1])
    rows = []
    for j in range(m+1):
        coeff = (f**j * modulus**(m-j)).coeffs()
        for shift in range(2):
            row = [0]*dim
            for i, value in enumerate(coeff, shift):
                row[i] = value*radius**i
            rows.append(row)
    vector = fmpz_mat(rows).lll().tolist()[0]
    poly = fmpz_poly([int(v)//radius**i for i, v in enumerate(vector)])
    for x, _ in poly.roots():
        j = center+int(x)
        if abs(int(x)) <= radius and (j*j-4*scale*scale) % modulus == 0:
            divisor = gcd(j-2*scale, n)
            if 1 < divisor < n:
                return divisor
    return None

def split_from_two_kloostermans(n, scale, first_sum, second_sum):
    """Raw sums at scale and H*scale, H=largest power of two ≤sqrt(n).

    Each supplied sum may have absolute error H/64.
    """
    from math import isqrt
    from mpmath import mp

    if n < 17 or n % 2 == 0 or gcd(scale, n) != 1:
        raise ValueError("Expected odd n >= 17 and a unit scale")
    modulus = n*n
    multiplier = 1 << (isqrt(n).bit_length()-1)
    with mp.workprec(4*n.bit_length()+64):
        angles = []
        for c, value in ((scale, first_sum), (multiplier*scale, second_sum)):
            cosine = mp.mpf(value)/(2*n)-mp.cospi(mp.mpf(4)*(c % modulus)/modulus)
            if abs(cosine) > mp.mpf("0.95"):
                return None
            angles.append(mp.acos(cosine)/(2*mp.pi))
        coarse, fine = angles
        centers = []
        for signed in (fine, -fine):
            winding = mp.nint(multiplier*coarse-signed)
            centers.append(int(mp.nint(modulus*(winding+signed)/multiplier)))
    for center in centers:
        divisor = _split_quadratic_interval(n, scale, center)
        if divisor is not None:
            return divisor
    return None


def split_from_radix_kloostermans(n, scale, evaluate, radix):
    """One trial using raw K_(n²)(1,c²) estimates within n/(128*radix)."""
    from mpmath import mp

    if n < 17 or n % 2 == 0 or gcd(scale, n) != 1:
        raise ValueError("Expected odd n >= 17 and a unit scale")
    if radix < 2 or radix & (radix-1):
        raise ValueError("Expected a power-of-two radix")
    modulus = n*n
    with mp.workprec(4*n.bit_length()+2*radix.bit_length()+64):
        def cosine(m):
            c = m*scale % modulus
            return mp.mpf(evaluate(c))/(2*n)-mp.cospi(mp.mpf(4)*c/modulus)

        ref = cosine(1)
        if abs(ref) > mp.mpf("0.75"):
            return None
        sine_ref = mp.sqrt(1-ref**2)
        angle = mp.acos(ref)/(2*mp.pi)
        power = 1
        while power < modulus:
            power *= radix
            direction = 1
            divisor = gcd(power+1, n)
            if 1 < divisor < n:
                return divisor
            if divisor == n:
                direction = -1
            value = cosine(power)
            sine = direction*(value*ref-cosine(power+direction))/sine_ref
            fraction = mp.arg(mp.mpc(value, sine))/(2*mp.pi)
            angle = (mp.nint(power*angle-fraction)+fraction)/power
        index = int(mp.nint(modulus*angle))
        divisor = gcd(index-2*scale, n)
        return divisor if 1 < divisor < n else None


def short_lift_mean_estimate(n, start, stride, length, precision=80):
    """Estimate a two-component mean; the error guarantee is over a random unit stride."""
    from mpmath import mp

    if n < 5 or gcd(n, 6*stride) != 1 or length < 1:
        raise ValueError("Expected a valid base, unit stride and positive length")
    with mp.workprec(precision+n.bit_length()):
        turn = 2j*mp.pi/n
        numerator = mp.expm1(turn*((stride*length) % n))
        denominator = length*mp.expm1(turn*(stride % n))
        phase = mp.expjpi(2*mp.mpf((1+n*start) % (n*n))/(n*n))
        return 4 + phase*numerator/denominator


def recover_from_fixed_length_windows(n, start, stride, length, window_mean):
    """Recover f(start) from length-L means on cyclic progressions modulo n."""
    if n < 1 or length < 1:
        raise ValueError("Expected a positive modulus and window length")
    start %= n
    result = window_mean(start, length*stride % n)
    for j in range(1, length):
        step = j*stride % n
        result += window_mean(start, step)-window_mean((start+step) % n, step)
    return result


def difference_from_window_means(n, left, right, length, window_mean, trials=9, draw=None):
    """Recover a real point difference from means, tolerating a minority of bad strides."""
    from secrets import randbelow
    from statistics import median

    if n < 2 or length < 1 or gcd(length, n) != 1 or trials < 1 or trials % 2 == 0:
        raise ValueError("Expected an invertible window length and a positive odd trial count")
    draw = randbelow if draw is None else draw
    left, right = left % n, right % n
    if left == right:
        return 0
    inverse = pow(length, -1, n)
    estimates = []
    for _ in range(trials):
        anchor = draw(n)
        def relative(start):
            step = (anchor-start)*inverse % n
            return window_mean(start, step)-window_mean((start+step) % n, step)
        estimates.append(length*(relative(left)-relative(right)))
    return median(estimates)
