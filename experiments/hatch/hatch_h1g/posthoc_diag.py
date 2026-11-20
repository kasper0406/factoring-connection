#!/usr/bin/env python3
"""Post-hoc diagnosis of P1c failure in E1c (does NOT modify pre-registered results.json)."""
import math, random, statistics
from toy_experiments import (
    primes_up_to, delta_p_iter, delta_p_mod, add_congruence_to_basis, hnf_reduce_rows, exact_lll
)

def diag():
    rng = random.Random(20261120)
    ps_odd = [q for q in primes_up_to(29) if q >= 3]
    P_odd = math.prod(ps_odd)
    print("ln(P_odd) =", round(math.log(P_odd), 4))

    for label, H_val, j_lo_frac in [
        ("Original H=2^10, j0 in [3, H)", 2**10, 0.0),
        ("Active s=3,4 regime: H=2^8 (b=5.55), P=prod_{3<=q<=79} q (W'=71.5), j0 in [3H/4, H)", 2**8, 0.75),
    ]:
        ps_use = ps_odd if H_val == 2**10 else [q for q in primes_up_to(79) if q >= 3]
        P_use = math.prod(ps_use)
        divs = [1]
        for q in ps_use:
            divs += [d * q for d in divs]
        divs = sorted(set(divs))
        # Use finer grid of 80 divisors around threshold
        cand = divs[::max(1, len(divs) // 80)]
        print(f"\n=== {label} (ln P={math.log(P_use):.2f}, b={math.log(H_val):.2f}) ===")
        for s in (2, 3, 4):
            buium_scales = [1] + [max(1, abs(delta_p_iter(H_val, 2, u))) for u in range(s)]
            copp_scales = [H_val**u for u in range(s + 1)]
            b_mins, c_mins, j0_vals = [], [], []
            for trial in range(20):
                lo = max(3, int(j_lo_frac * H_val))
                j0 = rng.randint(lo, H_val - 1)
                j0_vals.append(j0)
                gb = gc = math.log(P_use)
                for d in cand:
                    u = rng.randrange(1, P_use)
                    while math.gcd(u, P_use // d) != 1:
                        u = rng.randrange(1, P_use)
                    w = (j0 - d * u) % P_use
                    lnd = math.log(d)
                    if gb == math.log(P_use):
                        d_dim = s + 1
                        B0 = [[1 if i == j else 0 for j in range(d_dim)] for i in range(d_dim)]
                        coeffs = [1] + [delta_p_mod(w, 2, u_idx, P_use) for u_idx in range(s)]
                        B_lat, _ = add_congruence_to_basis(B0, coeffs, P_use)
                        B_lat = hnf_reduce_rows(B_lat)
                        B_sc = [[row[k] * buium_scales[k] for k in range(d_dim)] for row in B_lat]
                        red = exact_lll(B_sc)
                        c_vec = [red[0][k] // buium_scales[k] for k in range(d_dim)]
                        val = c_vec[0] + sum(c_vec[u_idx + 1] * delta_p_iter(j0, 2, u_idx) for u_idx in range(s))
                        if any(c_vec) and val == 0:
                            gb = lnd
                            if trial == 0:
                                print(f"  [s={s} Buium first hit] j0={j0}, lnd={lnd:.2f}, c_vec={c_vec}")
                    if gc == math.log(P_use):
                        d_dim = s + 1
                        B0 = [[1 if i == j else 0 for j in range(d_dim)] for i in range(d_dim)]
                        coeffs = [pow(w, u_idx, P_use) for u_idx in range(s + 1)]
                        B_lat, _ = add_congruence_to_basis(B0, coeffs, P_use)
                        B_lat = hnf_reduce_rows(B_lat)
                        B_sc = [[row[k] * copp_scales[k] for k in range(d_dim)] for row in B_lat]
                        red = exact_lll(B_sc)
                        c_vec = [red[0][k] // copp_scales[k] for k in range(d_dim)]
                        val = sum(c_vec[u_idx] * (j0**u_idx) for u_idx in range(s + 1))
                        if any(c_vec) and val == 0:
                            gc = lnd
                            if trial == 0:
                                print(f"  [s={s} Copp  first hit] j0={j0}, lnd={lnd:.2f}, c_vec={c_vec}")
                b_mins.append(gb)
                c_mins.append(gc)
            print(f"  s={s}: mean Buium={statistics.mean(b_mins):.3f}, mean Copp_r1={statistics.mean(c_mins):.3f}")

if __name__ == "__main__":
    diag()
