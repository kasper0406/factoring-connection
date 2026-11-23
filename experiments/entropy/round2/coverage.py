from inference import *


def main():
    times = list(csv.DictReader(open(RES / "timings.tsv"), delimiter="\t"))
    output = []
    for i in range(1, 8):
        name = f"G{i}"
        ts = [v for v in times if v["name"].startswith(name + "/")]
        eligible = [v for v in ts if v["eligible"] == "1"]
        output.append(
            dict(
                family=name,
                implemented=len(ts),
                eligible_at_64bits=len(eligible),
                frozen_unary_tests=sum(
                    1 for v in CAT[:NU] if v["expression"].startswith(name + "/")
                ),
                min_median_us=min(float(v["median_us"]) for v in eligible),
                max_median_us=max(float(v["median_us"]) for v in eligible),
            )
        )
    write_records("grammar_coverage.tsv", output)


if __name__ == "__main__":
    main()
