from screen import *
import sys
from concurrent.futures import ThreadPoolExecutor


def locks(n):
    rows = list(csv.DictReader(open(RES / "locks.tsv"), delimiter="\t"))
    return rows[: 5 if n >= 80 else 20]


def cfg(n):
    lines = []
    for d in locks(n):
        t = int(d["test"])
        op, a, b, c = map(int, OPS[t])
        un = [t] if op == 0 else ([a, b] if op == 1 else [a, b, c])
        vals = [str(op), str(len(un))]
        for i in un:
            vals.extend(map(str, OPS[i, 1:]))
        lines.append(" ".join(vals))
    (DATA / f"lock{n}.cfg").write_text("\n".join(lines) + "\n")


def jobs(n):
    if n <= 56:
        return [path(n, "C")]
    return sorted(DATA.glob(f"C{n}_*.pairs"))


def one_job(job):
    n, fn = job
    t = time.monotonic()
    tmp = Path(str(fn) + ".locked.tmp")
    dst = Path(str(fn) + ".locked")
    run(["./locked", n, fn.name, DATA / f"lock{n}.cfg", fn, tmp, "1"])
    assert tmp.stat().st_size == fn.stat().st_size
    tmp.replace(dst)
    print(n, fn.name, round(time.monotonic() - t, 2), flush=True)


if __name__ == "__main__":
    assert (R / "locked.utc").exists()
    for n in (40, 48, 56, 64, 80, 96):
        cfg(n)
    work = [(n, fn) for n in (40, 48, 56, 64, 96, 80) for fn in jobs(n)]
    with ThreadPoolExecutor(7) as ex:
        list(ex.map(one_job, work))
