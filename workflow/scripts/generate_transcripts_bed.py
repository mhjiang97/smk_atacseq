from pybedtools import BedTool

IDX_CHROM = 0
IDX_FEATURE = 2
IDX_START = 3
IDX_END = 4
IDX_STRAND = 6
IDX_ATTRIBUTES = 8


def parse_attributes(text):
    parsed = {}
    for item in text.strip().rstrip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        key, _, rest = item.partition(" ")
        parsed.setdefault(key, set()).add(rest.strip().strip('"'))
    return parsed


def matches_attributes(text, required):
    parsed = parse_attributes(text)
    return all(
        parsed.get(key, set()) & set(values) for key, values in required.items()
    )


def main():
    intervals = []
    with open(path_gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[IDX_FEATURE] != feature:
                continue
            if attributes and not matches_attributes(fields[IDX_ATTRIBUTES], attributes):
                continue
            start = int(fields[IDX_START]) - 1
            end = int(fields[IDX_END])
            if start >= end:
                continue
            intervals.append(
                f"{fields[IDX_CHROM]}\t{start}\t{end}\t.\t0\t{fields[IDX_STRAND]}"
            )

    if not intervals:
        raise ValueError(f"No {feature} records found in {path_gtf}.")

    merged = (
        BedTool("\n".join(intervals), from_string=True)
        .sort()
        .merge(s=True, c=6, o="distinct")
    )
    with open(path_bed, "w") as fh:
        fh.writelines(f"{iv.chrom}\t{iv.start}\t{iv.end}\t.\t0\t{iv.fields[-1]}\n" for iv in merged)


path_gtf = snakemake.input.gtf  # type: ignore # noqa: F821
path_bed = snakemake.output.bed  # type: ignore # noqa: F821
feature = snakemake.params.feature  # type: ignore # noqa: F821
attributes = snakemake.params.attributes  # type: ignore # noqa: F821

if __name__ == "__main__":
    main()
