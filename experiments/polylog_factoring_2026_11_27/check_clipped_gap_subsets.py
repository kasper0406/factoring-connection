"""Independent finite checks of the two-query clipped-gap SubsetSum encoding."""

from fractions import Fraction
from itertools import product
from math import gcd, prod

from sympy import nextprime

if __package__:
    from .clipped_gap_subsets import deterministic_components, encode, target_scale, upper_query
    from .check_single_quadratic_count import direct_contrast, packed_count
else:
    from clipped_gap_subsets import deterministic_components, encode, target_scale, upper_query
    from check_single_quadratic_count import direct_contrast, packed_count


def crt(residues, moduli):
    n = prod(moduli)
    return sum(x*(n//m)*pow(n//m,-1,m) for x,m in zip(residues,moduli)) % n


def stationary_table(instance):
    n, base = instance.modulus, instance.base
    squares = tuple(q*q for q in instance.components)
    arc = Fraction(2*(sum(instance.weights)+1)*n,instance.scale)
    assert 0 < instance.start < instance.start+instance.width < arc and 2*arc < n
    records = []
    for labels in product((1,-1), repeat=len(squares)):
        root = crt(labels,squares)
        lower = base*root % n
        assert min(lower,n-lower) < arc
        # CRT here does not use the suite's weighted-coefficient construction.
        assert lower == sum(e*a for e,a in zip(labels,instance.coefficients)) % n
        selected = instance.start <= lower < instance.start+instance.width
        weight = sum(w for w,e in zip(instance.weights,labels) if e == 1)
        assert selected == (weight == instance.target)
        assert ((lower+instance.offset) % n < instance.width) == selected
        if selected:
            assert labels[0] == 1 and (lower+instance.offset) % n <= (n-1)//2
            phases = []
            for digit in (0, *(q.upper_digit for q in instance.queries)):
                scale = base+n*digit
                phase = 2*scale*root % (n*n)
                # An independent lifted hyperbola point giving this phase.
                for y in (0,1,n-1):
                    x = lower+n*y
                    assert x*(phase-x) % (n*n) == scale**2 % (n*n)
                phases.append(phase)
            for phase,expected in zip(phases[1:],(Fraction(1,4),Fraction(3,4))):
                error = abs(Fraction(phase,n*n)-expected)
                assert error <= sum(Fraction(1,2*q) for q in instance.components)
            records.append(tuple(phase//n <= (n-1)//2 for phase in phases))
    return records


def off_root_mass(n, base, start, width, primes, exponents=None):
    """Exact counts by gcd class, even when the lower window is enormous."""
    exponents = (1,)*len(primes) if exponents is None else tuple(exponents)
    profiles = list(product(*(range(e+1) for e in exponents)))
    divisible = {}
    divisors = {}
    for profile in profiles:
        chosen = [p**e for p,e in zip(primes,profile) if e]
        divisor = prod(chosen)
        divisors[profile] = divisor
        count = 0
        for labels in product((1,-1),repeat=len(chosen)):
            a = base*crt(labels,chosen) % divisor
            count += (start+width-1-a)//divisor - (start-1-a)//divisor
        divisible[profile] = count
    # Successive differences invert divisibility counts on a product of chains.
    exact = divisible.copy()
    for j,e in enumerate(exponents):
        for v in profiles:
            if v[j] < e:
                upper = (*v[:j],v[j]+1,*v[j+1:])
                exact[v] -= exact[upper]
    assert min(exact.values()) >= 0 and sum(exact.values()) == width
    mass = sum(divisors[v]*count for v,count in exact.items() if 1 < divisors[v] < n)
    largest = max((divisors[v] for v,count in exact.items() if count and divisors[v] < n),
                  default=0)
    return mass,exact[exponents],largest


def check_case(weights, target, *, components=None, exact_mass=True):
    if components is None:
        k = len(weights)+1
        scale = target_scale(weights)
        first = 1000*k*scale
        primes = tuple(int(nextprime(first+1000*step+1)) for step in range(k))
        components = scale,primes,(1,)*k
    scale,primes,exponents = components
    instance = encode(weights,target,primes,exponents=exponents,scale=scale)
    records = stationary_table(instance)
    solutions = sum(sum(w for w,bit in zip(weights,mask) if bit) == target
                    for mask in product((0,1),repeat=len(weights)))
    assert len(records) == solutions
    reference = sum(row[0] for row in records)
    stationary = tuple(sum(row[i]-row[0] for row in records)*query.polarity
                       for i,query in enumerate(instance.queries,1))
    assert stationary == (solutions-reference,reference)
    n = instance.modulus
    if max(exponents) == 1:
        analytic = instance.width*3**len(primes) + prod(1+2*p for p in primes)-2**len(primes)*n
    else:
        # Near-zero stationary representatives exclude nonstationary divisors above this cap.
        largest_allowed = Fraction(8*(sum(instance.weights)+1)*n,scale)
        analytic = prod(1+2*e for e in exponents)*(instance.width+largest_allowed)
    assert 4*analytic < n
    if exact_mass:
        mass,count,largest = off_root_mass(n,instance.base,instance.start,instance.width,
                                          primes,exponents)
        assert count == solutions and mass <= analytic
        if max(exponents) > 1:
            assert largest <= largest_allowed
    else:
        mass = analytic
    # Mass bounds the entire unknown signed difference with unit lower weights.
    assert any(h*n-mass >= 3*Fraction(n,4) for h in stationary) == bool(solutions)
    if not solutions:
        assert all(h*n+mass <= Fraction(n,4) for h in stationary)
    ratio = Fraction(mass,n)
    shift = ratio.denominator.bit_length()-ratio.numerator.bit_length()-1
    bound_text = f"<2^-{shift}" if shift > 100 else f"{float(ratio):.3g}"
    print(f"weights={weights}, target={target}, {n.bit_length()}-bit N, "
          f"solutions={solutions}, reference={reference}, signed main terms={stationary}, "
          f"off-root bound/N={bound_text}, exact gcd histogram={exact_mass}")


def literal_controls():
    # Tiny original clipped queries: no promise on the size of the primes.
    # Compare complete direct hyperbola counts with the signed stationary contribution.
    n, base, primes = 15, 1, (3,5)
    width, start = 2, 3
    offset = -start % n
    raw_roots = [a for a in range(start,start+width) if (a*a-base*base) % n == 0]
    for q in (upper_query(base,primes,1),upper_query(base,primes,-1)):
        j = q.upper_digit
        coefficients = (n,base,j,width,offset,1,offset)
        contrast = direct_contrast(*coefficients)
        baseline = width*n*n*((n+1)//2)
        assert packed_count(*coefficients,q.polarity) == baseline+q.polarity*contrast
        stationary = 0
        for a in raw_roots:
            x = a
            phases = [(x+s*s*pow(x,-1,n*n)) % (n*n) for s in (base,base+n*j)]
            stationary += n*(int(phases[1]//n <= n//2)-int(phases[0]//n <= n//2))
        mass,count,_ = off_root_mass(n,base,start,width,primes)
        assert count == len(raw_roots) and abs(contrast-stationary) <= mass

    # This second toy has both off-root and stationary lower blocks.
    n,base,primes = 35,1,(5,7)
    start,width = 5,8
    mass,count,_ = off_root_mass(n,base,start,width,primes)
    assert mass == sum(g for a in range(start,start+width)
                       if 1 < (g := gcd(a*a-base*base,n)) < n) and count == 1
    seen = set()
    for j in range(n):
        offset = -start % n
        raw = direct_contrast(n,base,j,width,offset,1,offset)
        a = 6
        value = lambda c: ((a+c*c*pow(a,-1,n*n)) % (n*n))//n <= n//2
        stationary = n*(value(base+n*j)-value(base))
        assert abs(raw-stationary) <= mass
        seen.add(raw-stationary)
    assert any(seen)
    print("Two original packed comparisons and 35 further literal contrasts verified")


if __name__ == "__main__":
    for weights,target in [((2,4,6),5),((2,3,7),5),((1,2,3,4),5),
                           ((1,1,1,1,1),2)]:
        check_case(weights,target)
    for weights,target in [((2,4),3),((2,3),2),((1,1),1)]:
        check_case(weights,target,components=(1<<27,(3,5,7),(25,17,14)))
    weights = (1,2)
    small = deterministic_components(weights,clear_predivision=False)
    assert small[1] == (3,5,7)
    check_case(weights,2,components=small,exact_mass=False)
    safe = deterministic_components(weights)
    _,primes,exponents = safe
    bits = prod(p**e for p,e in zip(primes,exponents)).bit_length()
    trial_bound = (1 << 31)*(bits+1)
    assert all(p > trial_bound for p in primes)
    check_case(weights,2,components=safe,exact_mass=False)
    print(f"Predivision passed: smallest prime {min(primes)}, "
          f"trial bound {trial_bound}, {bits}-bit modulus")
    literal_controls()
