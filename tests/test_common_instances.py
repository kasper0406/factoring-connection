from math import prod

from src.common import instances as I
from src.common.baselines import pollard_rho_brent, trial_division
from src.common.metrics import Metrics, append_jsonl, read_jsonl


def test_balanced_bits_and_determinism():
    for b in I.BALANCED_BITS:
        a = I.balanced_semiprime(b, 3)
        assert a.N.bit_length() == b and prod(a.factors) == a.N
        assert all(I.is_prime(x) for x in a.factors) and a.factors[0] != a.factors[1]
        assert I.balanced_semiprime(b, 3).N == a.N
    assert I.balanced_semiprime(32, 0).N != I.balanced_semiprime(32, 1).N
    assert "labels" not in repr(I.balanced_semiprime(16, 0))


def test_other_families():
    for inst in [I.unbalanced_semiprime(64, 12, 0), I.prime_power(40, 3, 0),
                 I.multi_prime(48, 3, 0), I.multi_prime(64, 4, 1),
                 I.carmichael_chernick(8, 0), I.diag_near_equal(32, 0),
                 I.diag_smooth_pm1(48, 0)]:
        assert prod(inst.factors) == inst.N and all(I.is_prime(x) for x in inst.factors)
        if inst.family != "carmichael":
            assert inst.N.bit_length() == inst.params["bits"]
    assert I.is_prime(I.prime_instance(64, 0).N)


def test_carmichael_korselt():
    for i, N in enumerate(I.CARMICHAEL_SMALL):
        inst = I.carmichael_small(i)
        fs = inst.factors
        assert len(set(fs)) == len(fs) >= 3
        assert all((N - 1) % (p - 1) == 0 for p in fs)
    c = I.carmichael_chernick(10, 2)
    assert all((c.N - 1) % (p - 1) == 0 for p in c.factors)


def test_baselines(tmp_path):
    for b in (16, 24, 32, 40):
        inst = I.balanced_semiprime(b, 0)
        m = Metrics()
        assert trial_division(inst.N, m) == inst.factors[0]
        d = pollard_rho_brent(inst.N, 0, m)
        assert d in inst.factors
        rec = m.to_record(N=inst.N)
        append_jsonl(str(tmp_path / "x.jsonl"), rec)
    assert len(read_jsonl(str(tmp_path / "x.jsonl"))) == 4
    assert trial_division(10007) is None
