"""Fast (<1s) numerical audit of cadobench testbed parameters, leakage traps, and power calculations
for the guide_spec.md specification (no heavy compute)."""
import json
import math
import glob
import os

REPO = "/Volumes/git/factoring-connection"

def factor_small(n: int, limit: int = 2000000):
    n = abs(n)
    facs = {}
    d = 2
    while d * d <= n and d <= limit:
        while n % d == 0:
            facs[d] = facs.get(d, 0) + 1
            n //= d
        d = 3 if d == 2 else d + 2
    if n > 1:
        facs[n] = facs.get(n, 0) + 1
    return facs

print("=== 1. CADOBENCH TESTBED PARAMETER & LEAKAGE AUDIT ===")
metas = sorted(glob.glob(os.path.join(REPO, "experiments/raw/cadobench_testbed_c*/meta.json")),
               key=lambda p: int(p.split("_c")[-1].split("/")[0]))

audit = []
for mp in metas:
    with open(mp) as fh:
        m = json.load(fh)
    dig = m["digits"]
    f0 = [int(x) for x in m["f0_rational_low_to_high"]]
    f1 = [int(x) for x in m["f1_algebraic_low_to_high"]]
    cp = m["cado_params"]
    lpb0, lpb1 = int(cp["lpb0"]), int(cp["lpb1"])
    B0, B1 = 1 << lpb0, 1 << lpb1
    lim0, lim1 = int(cp["lim0"]), int(cp["lim1"])
    qmin = int(m["qmin"])
    g0, g1 = f0[0], f0[1]
    c0, cd = f1[0], f1[-1]
    fac_g1 = factor_small(g1, limit=B0)
    fac_cd = factor_small(cd, limit=B1)
    y = 2000
    g1_in_yB = {p: e for p, e in fac_g1.items() if y < p <= B0}
    cd_in_yB = {p: e for p, e in fac_cd.items() if y < p <= B1}
    print(f"c{dig}: lpb=({lpb0},{lpb1}) B=({B0},{B1}) lim=({lim0},{lim1}) qmin={qmin} (qmin<B1: {qmin < B1})")
    print(f"   Y1 (g1)={g1} factors={fac_g1} -> LEAKING PRIMES in (2000, B0]: {g1_in_yB}")
    print(f"   ad (cd)={cd} factors={fac_cd} -> LEAKING PRIMES in (2000, B1]: {cd_in_yB}")
    audit.append((dig, B0, B1, qmin, g1_in_yB, cd_in_yB))

print("\n=== 2. CUNNINGHAM PRIMES IN (2000, 2^21] (RAW-BIT LINEAR-FORM LEAKAGE) ===")
def is_prime(n):
    if n < 2: return False
    if n % 2 == 0: return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0: return False
        d += 2
    return True

cunningham = set()
for m in range(6, 25):
    for val in ((1 << m) - 1, (1 << m) + 1):
        for p in factor_small(val, 3000000):
            if 2000 < p <= (1 << 21) and is_prime(p):
                cunningham.add((p, m, "+" if val == (1 << m) + 1 else "-"))
cunn_sorted = sorted(cunningham)
print(f"Found {len(cunn_sorted)} prime divisors of 2^m +/- 1 (m<=24) in (2000, 2^21]:")
print("  First 15:", cunn_sorted[:15])

print("\n=== 3. PLANTED SIGNAL KL / LOG-LOSS / POWER CALIBRATION ===")
# Verify exact KL vs finite-p0 binary cross-entropy improvement via Gauss-Hermite quadrature
import math
# Standard normal nodes/weights via simple grid [-8, 8]
hs = 0.01
zs = [i * hs for i in range(-800, 801)]
phi = [math.exp(-0.5 * z * z) / math.sqrt(2 * math.pi) * hs for z in zs]

for u0, u1, p0_base in [(3.2, 4.4, 5e-4), (3.5, 4.8, 1e-4), (2.5, 3.2, 5e-3)]:
    ulogu = u0 * math.log(u0) + u1 * math.log(u1)
    print(f"\nOperating point (u0={u0}, u1={u1}), sum u*ln(u) = {ulogu:.3f}, base rate p0 = {p0_base:.1e}:")
    for theta in [0.01, 0.02, 0.05, 0.10, 0.25]:
        kl_target = theta * ulogu
        sigma = math.sqrt(2.0 * kl_target)
        # Marginal P(G=z) = (1-p0)*N(z;0,1) + p0*N(z;sigma,1)
        # Posterior p1(z) = p0*N(z;sigma,1) / P(G=z)
        ll_base = -p0_base * math.log(p0_base) - (1 - p0_base) * math.log(1 - p0_base)
        ll_mod = 0.0
        var_llr_hit = 0.0
        for z, w0 in zip(zs, phi):
            w1 = math.exp(-0.5 * (z - sigma) ** 2) / math.sqrt(2 * math.pi) * hs
            pg = (1 - p0_base) * w0 + p0_base * w1
            p1 = (p0_base * w1) / pg
            ll_mod += -(p0_base * w1 * math.log(p1) + (1 - p0_base) * w0 * math.log(1 - p1))
        dll = ll_base - ll_mod
        ihit = dll / p0_base
        # Under Y=1, z ~ N(sigma, 1), log(w1/w0) = sigma*z - 0.5*sigma^2 has mean 0.5*sigma^2 = kl_target and SD sigma!
        # With H smooth hits, SE(mean LLR on hits) approx sigma / sqrt(H); z-score = sqrt(H) * (0.5*sigma^2)/sigma = sqrt(H * 0.5 * kl_target)!
        # Wait: let's also include Y=0 variance: total per-candidate score has Fisher info / variance ~ 2 * p0 * kl_target!
        # Let's check: E[(dll_i)^2] where dll_i = y*log(p1/p0) + (1-y)*log((1-p1)/(1-p0)):
        m1 = sum((p0_base * (math.exp(-0.5*(z-sigma)**2)/math.sqrt(2*math.pi)*hs) * math.log(((p0_base * math.exp(-0.5*(z-sigma)**2)/math.sqrt(2*math.pi)*hs)/((1-p0_base)*w0 + p0_base*math.exp(-0.5*(z-sigma)**2)/math.sqrt(2*math.pi)*hs))/p0_base)
                  + (1-p0_base) * w0 * math.log((1 - (p0_base * math.exp(-0.5*(z-sigma)**2)/math.sqrt(2*math.pi)*hs)/((1-p0_base)*w0 + p0_base*math.exp(-0.5*(z-sigma)**2)/math.sqrt(2*math.pi)*hs))/(1-p0_base)))
                 for z, w0 in zip(zs, phi))
        m2 = sum((p0_base * (math.exp(-0.5*(z-sigma)**2)/math.sqrt(2*math.pi)*hs) * (math.log(((p0_base * math.exp(-0.5*(z-sigma)**2)/math.sqrt(2*math.pi)*hs)/((1-p0_base)*w0 + p0_base*math.exp(-0.5*(z-sigma)**2)/math.sqrt(2*math.pi)*hs))/p0_base))**2
                  + (1-p0_base) * w0 * (math.log((1 - (p0_base * math.exp(-0.5*(z-sigma)**2)/math.sqrt(2*math.pi)*hs)/((1-p0_base)*w0 + p0_base*math.exp(-0.5*(z-sigma)**2)/math.sqrt(2*math.pi)*hs))/(1-p0_base)))**2)
                 for z, w0 in zip(zs, phi))
        sd_i = math.sqrt(max(1e-30, m2 - m1 * m1))
        # For z_crit = 3.5 (Holm FWER 0.01 over ~20 tests) + z_power = 2.33 (99% power) -> total 5.83 SDs:
        n_req = (5.83 * sd_i / dll) ** 2
        hits_req = n_req * p0_base
        print(f"  theta={theta:4.2f}: target KL={kl_target:.4f} nats/hit, exact I_hit={ihit:.4f}, dLL={dll:.3e} nats/cand, H_99%pow(FWER.01)={hits_req:.0f} hits (N={n_req:.2e} cands)")
