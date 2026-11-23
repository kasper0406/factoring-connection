from confirm_eval import *
from budget import go


def worker(job):
    n, i = job
    fn = DATA / f"S{n}_{i:03d}.pairs"
    tmp = Path(str(fn) + ".partial")
    if not fn.exists():
        run(
            ["./sampler", n, 2, 1000000, f"power-addendum-{i}", tmp]
            + (["equal8"] if n == 96 else [])
        )
        tmp.replace(fn)
    go((n, fn))
    # At n=96 only the four locks defined on the sampled fiber are valid here.
    conf = DATA / f"lock{n}.cfg" if n == 64 else DATA / "sup96.cfg"
    out = Path(str(fn) + ".locked.tmp")
    run(["./locked", n, fn.name, conf, fn, out, "1"])
    assert out.stat().st_size == 16000000
    out.replace(str(fn) + ".locked")
    print(n, i, flush=True)


if __name__ == "__main__":
    assert (R / "power_addendum.utc").exists()
    (DATA / "sup96.cfg").write_text(
        "\n".join((DATA / "lock96.cfg").read_text().splitlines()[:4]) + "\n"
    )
    with ThreadPoolExecutor(6) as ex:
        list(
            ex.map(
                worker, [(64, i) for i in range(128)] + [(96, i) for i in range(170)]
            )
        )
