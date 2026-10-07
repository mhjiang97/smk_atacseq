def _get_samtools_view_args(wildcards):
    sample = wildcards.sample

    if sample in SAMPLES_PE:
        return "-f PAIRED -F UNMAP,MUNMAP,DUP,SECONDARY,SUPPLEMENTARY"
    if sample in SAMPLES_SE:
        return "-F UNMAP,DUP,SECONDARY,SUPPLEMENTARY"

    raise ValueError(f"Sample {sample} not found in SAMPLES_PE or SAMPLES_SE.")


rule samtools_view:
    input:
        bam="bowtie2/{sample}/{sample}.sorted.md.bam",
        bai="bowtie2/{sample}/{sample}.sorted.md.bam.bai",
        bed="ref/include.bed",
    output:
        bam=protected("bowtie2/{sample}/{sample}.sorted.md.filtered.bam"),
        bai=protected("bowtie2/{sample}/{sample}.sorted.md.filtered.bam.csi"),
    log:
        "logs/{sample}/samtools_view.log",
    threads: 1
    params:
        args=_get_samtools_view_args,
        name_mito="chrM",
        min_mapq=30,
    shell:
        """
        samtools view \\
            -@ {threads} \\
            -b \\
            -L {input.bed} \\
            -M \\
            {params.args} \\
            -q {params.min_mapq} \\
            -e 'rname != "{params.name_mito}"' \\
            --write-index \\
            -o {output.bam} \\
            {input.bam} \\
            1>{log} 2>&1
        """


def _get_samtools_view_split_args(wildcards):
    ranges = {
        "nfr": (0, 100),
        "1n": (180, 247),
        "2n": (315, 473),
        "3n": (558, 615),
    }
    type_fragment = wildcards.type_fragment
    if type_fragment not in ranges:
        raise ValueError(
            f"Type fragment {type_fragment} not found in nfr, 1n, 2n, 3n."
        )

    lo, hi = ranges[type_fragment]
    if lo == 0:
        return f"'(tlen >= -{hi} && tlen <= {hi}) && tlen != 0'"
    return f"'((tlen >= {lo} && tlen <= {hi}) || (tlen <= -{lo} && tlen >= -{hi}))'"


rule samtools_view_split:
    input:
        bam="bowtie2/{sample_pe}/{sample_pe}.sorted.md.filtered.bam",
    output:
        bam="bowtie2/{sample_pe}/{sample_pe}.sorted.md.filtered.{type_fragment}.bam",
    log:
        "logs/{sample_pe}/samtools_view_split.{type_fragment}.log",
    threads: 1
    params:
        args=_get_samtools_view_split_args,
    shell:
        """
        samtools view \\
            -b \\
            --write-index \\
            -e {params.args} \\
            -o {output.bam} \\
            {input.bam} \\
            1>{log} 2>&1
        """
