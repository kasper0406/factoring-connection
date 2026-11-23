"""Write input sizes and cryptographic digests in a replayable order."""

import csv
import hashlib
from config import DATA, R, RES


def main():
    samples = sorted([*DATA.glob("*.pairs"), *DATA.glob("*_[DVC].u64")])
    with (RES / "data_manifest.tsv").open("w") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_file", "draws", "bytes", "sha256"])
        for sample in samples:
            with sample.open("rb") as source:
                digest = hashlib.file_digest(source, "sha256").hexdigest()
            byte_count = sample.stat().st_size
            writer.writerow(
                [sample.relative_to(R), byte_count // 16, byte_count, digest]
            )


if __name__ == "__main__":
    main()
