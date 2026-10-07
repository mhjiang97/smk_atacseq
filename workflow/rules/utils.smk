def get_targets():
    targets = []

    targets += [
        f"bowtie2/{sample}/{sample}.sorted.md.filtered{ext}"
        for sample in SAMPLES
        for ext in [
            ".lc_extrap.txt",
            ".c_curve.txt",
            ".bound_pop.txt",
            ".quality_by_cycle.pdf",
        ]
    ]
    targets += [
        f"bowtie2/{sample}/{sample}.md.filtered.shifted.sorted.bw"
        for sample in SAMPLES_PE
    ]
    targets += [
        f"bedtools/{sample}/{sample}.md.filtered.shifted.extended.sorted.bw"
        for sample in SAMPLES_SE
    ]
    targets += [f"macs3/{sample}/{sample}_peaks.narrowPeak" for sample in SAMPLES]
    targets += [f"genrich/{sample}/{sample}_peaks.narrowPeak" for sample in SAMPLES]

    if not SAMPLES_PE.empty:
        targets += ["ataqv/mkarv.{caller}/index.html" for caller in CALLERS]

    targets += [f"encode/{sample}/{sample}.tss_enrich.qc" for sample in SAMPLES]
    targets += [
        f"frip/{sample}/{sample}.frip.{caller}.txt"
        for sample in SAMPLES
        for caller in CALLERS
    ]
    targets += [f"deeptools/{sample}/{sample}.heatmap.pdf" for sample in SAMPLES]
    targets += [
        f"deeptools/{sample}/{sample}.splitted.heatmap.pdf" for sample in SAMPLES_PE
    ]

    if TO_CLEAN_FQ:
        targets += [f"fastp/{sample}/{sample}.json" for sample in SAMPLES]
        if TO_RUN_FASTQC:
            targets += [f"fastqc/fastp/{sample}" for sample in SAMPLES]
            if TO_RUN_MULTIQC:
                targets += ["multiqc/fastp/multiqc_report.html"]
    else:
        if TO_RUN_FASTQC:
            targets += [f"fastqc/{sample}" for sample in SAMPLES]
            if TO_RUN_MULTIQC:
                targets += ["multiqc/multiqc_report.html"]

    return targets


def get_fastq_files(wildcards, named=True):
    sample = wildcards.sample
    dir_base = "fastp/{sample}" if TO_CLEAN_FQ else DIR_DATA

    if sample in SAMPLES_PE:
        files = {
            "fq1": f"{dir_base}/{{sample}}_R1.fq.gz",
            "fq2": f"{dir_base}/{{sample}}_R2.fq.gz",
        }
    elif sample in SAMPLES_SE:
        files = {"fq": f"{dir_base}/{{sample}}.fq.gz"}
    else:
        raise ValueError(f"Sample {sample} not found in SAMPLES_PE or SAMPLES_SE.")

    return files if named else list(files.values())


def get_extra_arguments(rule_name):
    if "args_extra" in config and rule_name in config["args_extra"]:
        return config["args_extra"][rule_name]

    return ""
