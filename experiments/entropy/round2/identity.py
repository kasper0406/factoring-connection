from confirm_eval import *


def audit_size(n):
    files0 = jobs(n)
    extra = (
        [DATA / f"S{n}_{i:03d}.pairs" for i in range(128 if n == 64 else 170)]
        if n in (64, 96)
        else []
    )
    while any(not Path(str(f) + ".locked").exists() for f in extra):
        time.sleep(3)
    manifest = DATA / f"identities{n}.list"
    manifest.write_text(
        "".join(f"{e} {f}\n" for e, fs in [(0, files0), (1, extra)] for f in fs)
    )
    with open(RES / f"identities{n}.json", "w") as out:
        run(["./identity", n, manifest, DATA / f"id{n}"], stdout=out)
    print(n, "identity audit", flush=True)


if __name__ == "__main__":
    for n in (40, 48, 56, 80, 64, 96):
        audit_size(n)
