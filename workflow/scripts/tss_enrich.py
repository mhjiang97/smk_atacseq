"""ENCODE-style TSS enrichment for ATAC-seq BAM files.

Greenleaf / ENCODE method: pile shifted reads in a +/- 2000 bp window around
each TSS, normalize by the mean of the outer 100 bp flanks, and take the
maximum of the profile as the enrichment score.

Works on single-end and paired-end BAMs (uses read coverage, not insert size).
"""

import logging
import os
from concurrent.futures import ProcessPoolExecutor

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pysam

log = logging.getLogger(__name__)

BP_EDGE = 2000
BINS = 400
FLANK_BP = 100


def load_chrom_sizes(path):
    sizes = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            sizes[fields[0]] = int(fields[1])
    return sizes


def load_tss(path):
    sites = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            start = int(fields[1])
            strand = fields[5] if len(fields) > 5 else "+"
            sites.append((chrom, start, strand))
    return sites


def detect_read_len(bam_path, n=10000):
    lengths = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            if read.is_unmapped or read.query_length is None:
                continue
            lengths.append(read.query_length)
            if len(lengths) >= n:
                break
    if not lengths:
        raise ValueError(f"Could not determine read length from {bam_path}.")
    return int(np.median(lengths))


def ensure_bam_index(bam_path):
    if os.path.exists(bam_path + ".bai") or os.path.exists(bam_path + ".csi"):
        return
    pysam.index(bam_path)


def _profile_chunk(args):
    bam_path, sites, chrom_sizes, bins, bp_edge, shift = args
    bin_width = (2 * bp_edge) / bins
    profiles = np.zeros((len(sites), bins), dtype=np.float64)
    bam = pysam.AlignmentFile(bam_path, "rb")
    refs = set(bam.references)

    for i, (chrom, tss, strand) in enumerate(sites):
        if chrom not in refs:
            continue
        chrom_len = chrom_sizes.get(chrom)
        window_start = tss - bp_edge
        window_end = tss + bp_edge
        fetch_start = max(0, window_start + min(shift, 0))
        fetch_end = window_end + max(shift, 0)
        if chrom_len is not None:
            fetch_end = min(fetch_end, chrom_len)

        for read in bam.fetch(chrom, fetch_start, max(fetch_start + 1, fetch_end)):
            if (
                read.is_unmapped
                or read.is_secondary
                or read.is_supplementary
                or read.is_duplicate
                or read.is_qcfail
            ):
                continue
            if strand == "+" and read.is_reverse:
                continue
            if strand == "-" and not read.is_reverse:
                continue

            rstart = read.reference_start + shift
            rend = read.reference_end + shift
            overlap_start = max(rstart, window_start)
            overlap_end = min(rend, window_end)
            if overlap_start >= overlap_end:
                continue
            b0 = int((overlap_start - window_start) / bin_width)
            b1 = int(np.ceil((overlap_end - window_start) / bin_width))
            b0 = max(0, b0)
            b1 = min(bins, b1)
            if b0 < b1:
                profiles[i, b0:b1] += 1

        if strand == "-":
            profiles[i] = profiles[i, ::-1]

    bam.close()
    return profiles


def _chunks(items, n):
    n = max(1, min(n, len(items)))
    size, extra = divmod(len(items), n)
    out = []
    start = 0
    for i in range(n):
        end = start + size + (1 if i < extra else 0)
        if start < end:
            out.append(items[start:end])
        start = end
    return out


def tss_array(bam_path, sites, chrom_sizes, bins, bp_edge, shift, threads):
    if not sites:
        raise ValueError("No TSS records found.")
    threads = max(1, threads)
    if threads == 1:
        return _profile_chunk((bam_path, sites, chrom_sizes, bins, bp_edge, shift))

    payload = [
        (bam_path, chunk, chrom_sizes, bins, bp_edge, shift)
        for chunk in _chunks(sites, threads)
    ]
    with ProcessPoolExecutor(max_workers=len(payload)) as pool:
        parts = list(pool.map(_profile_chunk, payload))
    return np.vstack(parts)


def make_tss_plot(bam_array, out_pdf, out_pdf_large, out_qc, bins=BINS, bp_edge=BP_EDGE):
    num_edge_bins = int(FLANK_BP / (2 * bp_edge / bins))
    bin_means = bam_array.mean(axis=0)
    avg_noise = (
        bin_means[:num_edge_bins].sum() + bin_means[-num_edge_bins:].sum()
    ) / (2 * num_edge_bins)
    if avg_noise == 0:
        raise ValueError("TSS flank noise is 0; check BAM chromosome names vs TSS BED.")
    bam_array = bam_array / avg_noise
    profile = bam_array.mean(axis=0)
    tss_point_val = float(profile.max())

    os.makedirs(os.path.dirname(out_qc) or ".", exist_ok=True)
    with open(out_qc, "w") as fh:
        fh.write(f"{tss_point_val}\n")

    x = np.linspace(-bp_edge, bp_edge, bins)
    fig, ax = plt.subplots()
    ax.plot(x, profile, color="r", label="Mean")
    ax.axvline(0, linestyle=":", color="k")
    ax.set_xlabel("Distance from TSS (bp)")
    ax.set_ylabel("TSS Enrichment")
    ax.legend(loc="best")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

    upper_prct = 99
    flat = bam_array.ravel()
    if np.percentile(flat, upper_prct) == 0.0:
        upper_prct = 100.0
    vmin = np.percentile(flat, 5)
    vmax = np.percentile(flat, upper_prct)
    order = np.argsort(bam_array.mean(axis=1))
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(5, 10),
        gridspec_kw={"height_ratios": [1, 4]},
        sharex=True,
    )
    axes[0].plot(x, profile, color="k")
    axes[0].fill_between(x, profile, color="k", alpha=0.3)
    axes[0].axvline(0, linestyle=":", color="k")
    axes[0].set_ylabel("TSS Enrichment")
    im = axes[1].imshow(
        bam_array[order],
        aspect="auto",
        interpolation="nearest",
        extent=[-bp_edge, bp_edge, 0, bam_array.shape[0]],
        vmin=vmin,
        vmax=vmax,
        cmap="viridis",
    )
    axes[1].set_xlabel("Distance from TSS (bp)")
    axes[1].set_ylabel("TSS (sorted by mean signal)")
    fig.colorbar(im, ax=axes[1], fraction=0.046, pad=0.04)
    fig.savefig(out_pdf_large, bbox_inches="tight")
    plt.close(fig)

    return tss_point_val


def main():
    logging.basicConfig(
        format="[%(asctime)s %(levelname)s] %(message)s",
        filename=path_log,
        filemode="w",
        level=logging.INFO,
    )
    ensure_bam_index(bam)
    chrom_sizes = load_chrom_sizes(chrsz)
    sites = load_tss(tss)
    length = read_len if read_len is not None else detect_read_len(bam)
    shift = -length / 2
    log.info(
        "TSS sites=%d read_len=%d shift=%.1f threads=%d",
        len(sites),
        length,
        shift,
        threads,
    )

    bam_array = tss_array(
        bam,
        sites,
        chrom_sizes,
        BINS,
        BP_EDGE,
        shift,
        threads,
    )
    score = make_tss_plot(bam_array, out_pdf, out_pdf_large, out_qc)
    log.info("TSS enrichment=%.4f", score)


bam = snakemake.input.bam  # type: ignore # noqa: F821
tss = snakemake.input.tss  # type: ignore # noqa: F821
chrsz = snakemake.input.fai  # type: ignore # noqa: F821
out_pdf = snakemake.output.pdf  # type: ignore # noqa: F821
out_pdf_large = snakemake.output.pdf_large  # type: ignore # noqa: F821
out_qc = snakemake.output.qc  # type: ignore # noqa: F821
path_log = snakemake.log[0]  # type: ignore # noqa: F821
threads = snakemake.threads  # type: ignore # noqa: F821
read_len = snakemake.params.len_read  # type: ignore # noqa: F821

if __name__ == "__main__":
    main()
