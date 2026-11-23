#!/usr/bin/env python3
"""Quantitative verification and adversarial audit script for entropy_spec.md.
Runs fast (< 5s), read-only wrt repo, self-contained in codex/scratch/ent_g/.
"""
from __future__ import annotations

import json
import math
import random
from statistics import NormalDist

ND = NormalDist()


def sieve_primes(lo: int, hi: int) -> list[int]:
    """Primes in [lo, hi)."""
    if hi <= 2:
        return []
    sieve = bytearray([1]) * hi
    sieve[0] = sieve[1] = 0
    for i in range(2, int(hi**0.5) + 1):
        if sieve[i]:
            sieve[i * i : hi : i] = bytearray(len(sieve[i * i : hi : i]))
    return [x for x in range(max(2, lo), hi) if sieve[x]]


def pnt_li_interval(lo: float, hi: float, steps: int = 200) -> float:
    """Simpson approximation to integral_{lo}^{hi} dx / ln(x)."""
    h = (hi - lo) / steps
    s = 1.0 / math.log(lo) + 1.0 / math.log(hi)
    for i in range(1, steps):
        x = lo + i * h
        s += (4.0 if i % 2 == 1 else 2.0) / math.log(x)
    return s * h / 3.0


def audit_universe_and_reuse() -> dict:
    sizes = [24, 32, 40, 48, 56, 64, 80, 96]
    m_total = 600_000  # 2e5 x 3 splits (disc/val/conf)
    out = {}
    for n in sizes:
        b = n // 2
        lo_top2, hi = 3 << (b - 2), 1 << b
        lo_full = 1 << (b - 1)
        if b <= 16:
            p_top2 = len(sieve_primes(lo_top2, hi))
            p_full = len(sieve_primes(lo_full, hi))
        else:
            p_top2 = round(pnt_li_interval(float(lo_top2), float(hi)))
            p_full = round(pnt_li_interval(float(lo_full), float(hi)))
        pairs_top2 = p_top2 * (p_top2 - 1) // 2
        pairs_full = p_full * (p_full - 1) // 2
        reuse_top2 = (2.0 * m_total) / p_top2
        out[n] = {
            "b": b,
            "primes_top2": p_top2,
            "primes_full": p_full,
            "pairs_top2": pairs_top2,
            "pairs_full": pairs_full,
            "prime_reuse_in_6e5_top2": round(reuse_top2, 3),
            "approx_unique_fraction_6e5_top2": round(
                (pairs_top2 / m_total) * (1.0 - math.exp(-m_total / pairs_top2))
                if m_total / pairs_top2 < 50
                else pairs_top2 / m_total,
                4,
            ),
        }
    return out


def audit_targets_and_biases() -> dict:
    """Exact census at b=12 (n=24) and sampled/exact checks at b=16 (n=32)."""
    res = {}
    for label, b, lo, hi in [
        ("n24_full", 12, 1 << 11, 1 << 12),
        ("n24_top2", 12, 3 << 10, 1 << 12),
    ]:
        ps = sieve_primes(lo, hi)
        pairs = [(ps[i], ps[j]) for i in range(len(ps)) for j in range(i + 1, len(ps))]
        m = len(pairs)
        bit_means = {}
        for off in range(1, 7):
            bit_pos = b - 1 - off
            if bit_pos >= 1:
                bit_means[f"top_minus_{off}"] = round(
                    sum((p >> bit_pos) & 1 for p, _ in pairs) / m, 5
                )
        n1m4 = [(p, q) for p, q in pairs if (p * q) % 4 == 1]
        eps_rate = sum(1 for p, _ in n1m4 if p % 4 == 1) / len(n1m4)
        n1m8 = [(p, q) for p, q in pairs if (p * q) % 8 == 1]
        qr2_rate = sum(1 for p, _ in n1m8 if p % 8 in (1, 7)) / len(n1m8)

        lin_1, loglin_1, pit_1 = [], [], []
        for p, q in pairs:
            N = p * q
            pmin = max(float(lo), N / float(hi))
            sqN = math.sqrt(N)
            if sqN - pmin <= 1e-9:
                continue
            u_lin = (p - pmin) / (sqN - pmin)
            u_loglin = 1.0 - math.log(q / p) / math.log(N / (pmin * pmin))
            smin = math.log(pmin)
            lnN = math.log(N)
            u_pit = 1.0 - math.log(math.log(q) / math.log(p)) / math.log((lnN - smin) / smin)
            lin_1.append(1 if u_lin >= 0.5 else 0)
            loglin_1.append(1 if u_loglin >= 0.5 else 0)
            pit_1.append(1 if u_pit >= 0.5 else 0)

        res[label] = {
            "n_primes": len(ps),
            "n_pairs": m,
            "sorted_bit_means": bit_means,
            "eps_pos_rate": round(eps_rate, 5),
            "qr2_pos_rate": round(qr2_rate, 5),
            "mean_u_lin_ge_half": round(sum(lin_1) / len(lin_1), 5),
            "mean_u_loglin_ge_half": round(sum(loglin_1) / len(loglin_1), 5),
            "mean_u_pit_ge_half": round(sum(pit_1) / len(pit_1), 5),
        }

    loglin_bias_by_n = {}
    for n in [24, 32, 48, 64, 96]:
        S = n * math.log(2.0)
        d = 0.5 * math.log(2.0)
        p_true = 1.0 - math.log((S / 2 + d / 2) / (S / 2 - d / 2)) / math.log(
            (S / 2 + d) / (S / 2 - d)
        )
        s_lin_mid = S / 2 - math.log(2.0 * math.sqrt(2.0) / (1.0 + math.sqrt(2.0)))
        p_lin_true = 1.0 - math.log((S - s_lin_mid) / s_lin_mid) / math.log(
            (S / 2 + d) / (S / 2 - d)
        )
        loglin_bias_by_n[n] = {
            "linear_adv_at_midN": round(abs(p_lin_true - 0.5), 6),
            "loglinear_adv_at_midN": round(abs(p_true - 0.5), 7),
            "one_over_n2": round(1.0 / (n * n), 7),
        }
    res["analytic_midpoint_bias"] = loglin_bias_by_n

    ps16 = sieve_primes(49000, 65600)
    in_range = [p for p in ps16 if 49152 <= p < 65536]
    unif_1m3 = sum(1 for p in in_range if p % 3 == 1) / len(in_range)
    unif_smooth = sum(1 for p in in_range if (p - 1) % 6 == 0) / len(in_range)
    prev_map = {ps16[i]: ps16[i - 1] for i in range(1, len(ps16))}
    tot_w = sum(p - prev_map[p] for p in in_range)
    gap_1m3 = sum((p - prev_map[p]) for p in in_range if p % 3 == 1) / tot_w
    gap_smooth = sum((p - prev_map[p]) for p in in_range if (p - 1) % 6 == 0) / tot_w

    def smooth13(x: int) -> int:
        s = 1
        for r in (2, 3, 5, 7, 11, 13):
            while x % r == 0:
                s *= r
                x //= r
        return s

    unif_s13 = sum(math.log2(smooth13(p - 1)) for p in in_range) / len(in_range)
    gap_s13 = sum((p - prev_map[p]) * math.log2(smooth13(p - 1)) for p in in_range) / tot_w
    res["nextprime_gap_bias_b16"] = {
        "n_primes": len(in_range),
        "unif_p1mod3": round(unif_1m3, 5),
        "nextprime_p1mod3": round(gap_1m3, 5),
        "unif_p1mod6": round(unif_smooth, 5),
        "nextprime_p1mod6": round(gap_smooth, 5),
        "unif_mean_log2_smooth13_pm1": round(unif_s13, 4),
        "nextprime_mean_log2_smooth13_pm1": round(gap_s13, 4),
    }
    return res


def audit_sample_sizes_and_power() -> dict:
    """Exact sample sizes for advantage 1/n^2 and 1/n, especially at n=96."""
    z_pow80 = ND.inv_cdf(0.80)
    z_pow90 = ND.inv_cdf(0.90)
    configs = {
        "single_1sided_a0.01_pow80": (ND.inv_cdf(0.99), z_pow80),
        "single_1sided_a0.01_pow90": (ND.inv_cdf(0.99), z_pow90),
        "single_2sided_a0.05_pow80": (ND.inv_cdf(0.975), z_pow80),
        "survivors20_1sided_Holm0.01_pow80": (ND.inv_cdf(1.0 - 0.01 / 20), z_pow80),
        "discovery_2e5_2sided_FWER0.01_pow80": (ND.inv_cdf(1.0 - 0.01 / (2 * 200_000)), z_pow80),
        "discovery_1e6_2sided_FWER0.01_pow80": (ND.inv_cdf(1.0 - 0.01 / (2 * 1_000_000)), z_pow80),
    }
    n96 = {}
    delta_n2 = 1.0 / (96 * 96)
    delta_n1 = 1.0 / 96.0
    for name, (zc, zp) in configs.items():
        m_n2 = math.ceil(((zc + zp) ** 2) / (4.0 * delta_n2 * delta_n2))
        m_n1 = math.ceil(((zc + zp) ** 2) / (4.0 * delta_n1 * delta_n1))
        n96[name] = {
            "z_crit": round(zc, 4),
            "z_sum": round(zc + zp, 4),
            "m_for_1_over_n2": m_n2,
            "m_for_1_over_n": m_n1,
        }

    m_codex = 200_000
    codex_2e5 = {}
    zc_1s = ND.inv_cdf(0.99)
    for n in [24, 32, 40, 48, 56, 64, 80, 96]:
        se = 0.5 / math.sqrt(m_codex)
        ez_n2 = (1.0 / (n * n)) / se
        ez_n1 = (1.0 / n) / se
        pow_n2 = 1.0 - ND.cdf(zc_1s - ez_n2)
        pow_n1 = 1.0 - ND.cdf(zc_1s - ez_n1)
        codex_2e5[n] = {
            "SE_adv": round(se, 6),
            "E_z_1_over_n2": round(ez_n2, 3),
            "power_1s_a0.01_1_over_n2": round(pow_n2, 4),
            "E_z_1_over_n": round(ez_n1, 2),
            "power_1s_a0.01_1_over_n": round(pow_n1, 4),
            "m_needed_single_1s_a0.01_pow80_1_over_n2": math.ceil(
                ((zc_1s + z_pow80) ** 2) * (n**4) / 4.0
            ),
        }

    mi_true_nats = 2.0 * (delta_n2**2)
    mi_null_bias_2e5_nats = 1.0 / (2.0 * m_codex)
    mi_null_bias_2e8_nats = 1.0 / (2.0 * 213_100_000)

    return {
        "n96_sample_sizes": n96,
        "codex_m2e5_power_table": codex_2e5,
        "mi_audit_n96": {
            "delta": round(delta_n2, 8),
            "true_MI_nats": mi_true_nats,
            "null_plugin_MI_bias_m2e5_nats": mi_null_bias_2e5_nats,
            "bias_to_signal_ratio_m2e5": round(mi_null_bias_2e5_nats / mi_true_nats, 1),
            "null_plugin_MI_bias_m2.13e8_nats": mi_null_bias_2e8_nats,
        },
    }


def fit_poly_wls(ns: list[int], y: list[float], w: list[float]) -> tuple[float, float, float]:
    """Grid profile WLS for C * n^(-k), C >= 0, k in [0.2, 6.0]."""
    best = (1e30, 0.0, 1.0)
    for gi in range(291):
        k = 0.2 + gi * 0.02
        basis = [n ** (-k) for n in ns]
        denom = sum(wi * bi * bi for wi, bi in zip(w, basis))
        num = sum(wi * yi * bi for wi, yi, bi in zip(w, y, basis))
        c = max(0.0, num / denom)
        chi2 = sum(wi * (yi - c * bi) ** 2 for wi, yi, bi in zip(w, y, basis))
        if chi2 < best[0]:
            best = (chi2, c, k)
    return best


def fit_exp_wls(ns: list[int], y: list[float], w: list[float]) -> tuple[float, float, float]:
    """Grid profile WLS for A * 2^(-beta * n), A >= 0, beta in [0.02, 0.60]."""
    best = (1e30, 0.0, 0.1)
    for gi in range(291):
        beta = 0.02 + gi * 0.002
        basis = [2.0 ** (-beta * n) for n in ns]
        denom = sum(wi * bi * bi for wi, bi in zip(w, basis))
        num = sum(wi * yi * bi for wi, yi, bi in zip(w, y, basis))
        a = max(0.0, num / denom)
        chi2 = sum(wi * (yi - a * bi) ** 2 for wi, yi, bi in zip(w, y, basis))
        if chi2 < best[0]:
            best = (chi2, a, beta)
    return best


def fit_exp_barrier_wls(ns: list[int], y: list[float], w: list[float]) -> tuple[float, float, float]:
    """Profile WLS for known barrier exponentials A * 2^(-beta * n), beta >= 1/6 (~0.1667)."""
    best = (1e30, 0.0, 1.0 / 6.0)
    for gi in range(251):
        beta = (1.0 / 6.0) + gi * 0.002
        basis = [2.0 ** (-beta * n) for n in ns]
        denom = sum(wi * bi * bi for wi, bi in zip(w, basis))
        num = sum(wi * yi * bi for wi, yi, bi in zip(w, y, basis))
        a = max(0.0, num / denom)
        chi2 = sum(wi * (yi - a * bi) ** 2 for wi, yi, bi in zip(w, y, basis))
        if chi2 < best[0]:
            best = (chi2, a, beta)
    return best


def fit_sub_wls(ns: list[int], y: list[float], w: list[float]) -> tuple[float, float, float]:
    """Grid profile WLS for A * exp(-gamma * sqrt(n * ln n)), A >= 0, gamma in [0.05, 2.5]."""
    best = (1e30, 0.0, 0.5)
    for gi in range(246):
        gamma = 0.05 + gi * 0.01
        basis = [math.exp(-gamma * math.sqrt(n * math.log(n))) for n in ns]
        denom = sum(wi * bi * bi for wi, bi in zip(w, basis))
        num = sum(wi * yi * bi for wi, yi, bi in zip(w, y, basis))
        a = max(0.0, num / denom)
        chi2 = sum(wi * (yi - a * bi) ** 2 for wi, yi, bi in zip(w, y, basis))
        if chi2 < best[0]:
            best = (chi2, a, gamma)
    return best


def audit_decay_fitting_protocol() -> dict:
    """Demonstrate the log-|adv| OLS trap vs our signed WLS + large-n gate + exact curvature SNR."""
    rng = random.Random(20261123)
    ns = [32, 40, 48, 56, 64, 80, 96]

    w_unit1e8 = [4.0 * 100_000_000.0] * len(ns)
    y_true_n2 = [1.0 / (n * n) for n in ns]
    ncp_barrier_per_1e8, _, beta_b = fit_exp_barrier_wls(ns, y_true_n2, w_unit1e8)
    ncp_free_exp_per_1e8, _, beta_f = fit_exp_wls(ns, y_true_n2, w_unit1e8)
    ncp_sub_per_1e8, _, gamma_s = fit_sub_wls(ns, y_true_n2, w_unit1e8)

    m_small = [200_000] * len(ns)
    ses_small = [0.5 / math.sqrt(m) for m in m_small]
    w_small = [1.0 / (se * se) for se in ses_small]

    naive_log_ols_poly_wins = 0
    naive_k_estimates = []
    our_protocol_fp = 0
    trials = 500
    x_logn = [math.log(n) for n in ns]
    mx = sum(x_logn) / len(ns)
    sxx_ln = sum((x - mx) ** 2 for x in x_logn)
    mn = sum(ns) / len(ns)
    sxx_n = sum((n - mn) ** 2 for n in ns)

    for _ in range(trials):
        y = [0.05 * (2.0 ** (-(n - 32) / 4.0)) + rng.gauss(0.0, se) for n, se in zip(ns, ses_small)]
        ly = [math.log(max(abs(yi), 1e-12)) for yi in y]
        my = sum(ly) / len(ly)
        b_poly = sum((x - mx) * (yi - my) for x, yi in zip(x_logn, ly)) / sxx_ln
        rss_poly = sum((yi - (my + b_poly * (x - mx))) ** 2 for x, yi in zip(x_logn, ly))
        b_exp = sum((n - mn) * (yi - my) for n, yi in zip(ns, ly)) / sxx_n
        rss_exp = sum((yi - (my + b_exp * (n - mn))) ** 2 for n, yi in zip(ns, ly))
        if rss_poly < rss_exp:
            naive_log_ols_poly_wins += 1
        naive_k_estimates.append(-b_poly)

        z80 = y[-2] / ses_small[-2]
        z96 = y[-1] / ses_small[-1]
        if z80 >= 3.0 and z96 >= 3.0:
            chi2_p, _, _ = fit_poly_wls(ns, y, w_small)
            chi2_eb, _, _ = fit_exp_barrier_wls(ns, y, w_small)
            if chi2_eb - chi2_p >= 9.21 and chi2_p <= 15.09:
                our_protocol_fp += 1

    m_powered = [100_000_000, 100_000_000, 100_000_000, 100_000_000, 100_000_000, 320_000_000, 320_000_000]
    ses_pow = [0.5 / math.sqrt(m) for m in m_powered]
    w_pow = [1.0 / (se * se) for se in ses_pow]

    poly_tp_vs_barrier = 0
    poly_tp_vs_free_exp = 0
    barrier_fp_powered = 0
    for _ in range(trials):
        y_p = [1.0 / (n * n) + rng.gauss(0.0, se) for n, se in zip(ns, ses_pow)]
        if y_p[-2] / ses_pow[-2] >= 3.0 and y_p[-1] / ses_pow[-1] >= 3.0:
            chi2_p, _, _ = fit_poly_wls(ns, y_p, w_pow)
            chi2_eb, _, _ = fit_exp_barrier_wls(ns, y_p, w_pow)
            chi2_ef, _, _ = fit_exp_wls(ns, y_p, w_pow)
            if chi2_eb - chi2_p >= 9.21 and chi2_p <= 15.09:
                poly_tp_vs_barrier += 1
            if chi2_ef - chi2_p >= 6.63 and chi2_p <= 15.09:
                poly_tp_vs_free_exp += 1

        y_b = [0.05 * (2.0 ** (-(n - 32) / 6.0)) + rng.gauss(0.0, se) for n, se in zip(ns, ses_pow)]
        if y_b[-2] / ses_pow[-2] >= 3.0 and y_b[-1] / ses_pow[-1] >= 3.0:
            chi2_p, _, _ = fit_poly_wls(ns, y_b, w_pow)
            chi2_eb, _, _ = fit_exp_barrier_wls(ns, y_b, w_pow)
            if chi2_eb - chi2_p >= 9.21 and chi2_p <= 15.09:
                barrier_fp_powered += 1

    return {
        "noiseless_ncp_for_1_over_n2_at_m1e8_per_size": {
            "vs_barrier_exp_beta_ge_1_6": round(ncp_barrier_per_1e8, 2),
            "closest_barrier_beta": round(beta_b, 4),
            "vs_unconstrained_exp": round(ncp_free_exp_per_1e8, 2),
            "closest_free_beta": round(beta_f, 4),
            "vs_L_half_subexp": round(ncp_sub_per_1e8, 2),
            "closest_subexp_gamma": round(gamma_s, 4),
        },
        "naive_log_abs_ols_on_2e5_N_minus_quarter_leak": {
            "poly_beats_exp_fraction": round(naive_log_ols_poly_wins / trials, 3),
            "mean_fitted_k": round(sum(naive_k_estimates) / len(naive_k_estimates), 2),
        },
        "signed_wls_protocol": {
            "fp_rate_on_N_minus_quarter_m2e5": round(our_protocol_fp / trials, 4),
            "fp_rate_on_N_minus_sixth_powered": round(barrier_fp_powered / trials, 4),
            "tp_rate_1_over_n2_vs_barrier_exp_powered": round(poly_tp_vs_barrier / trials, 3),
            "tp_rate_1_over_n2_vs_free_exp_powered": round(poly_tp_vs_free_exp / trials, 3),
        },
    }


def main() -> None:
    report = {
        "universe_and_reuse": audit_universe_and_reuse(),
        "targets_and_biases": audit_targets_and_biases(),
        "sample_sizes_and_power": audit_sample_sizes_and_power(),
        "decay_fitting_protocol": audit_decay_fitting_protocol(),
    }
    out_path = "/Volumes/git/factoring-connection/codex/scratch/ent_g/spec_audit_numbers.json"
    with open(out_path, "w") as fh:
        json.dump(report, fh, indent=2)
    print(json.dumps(report["decay_fitting_protocol"], indent=2))


if __name__ == "__main__":
    main()
