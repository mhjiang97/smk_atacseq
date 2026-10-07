files_bowtie2_index = multiext(
    config["index_bowtie2"],
    ".1.bt2",
    ".2.bt2",
    ".3.bt2",
    ".4.bt2",
    ".rev.1.bt2",
    ".rev.2.bt2",
)


def _get_fastq_input_arg(wildcards):
    sample = wildcards.sample
    dir_base = f"fastp/{sample}" if TO_CLEAN_FQ else DIR_DATA
    if sample in SAMPLES_PE:
        arg = f"-1 {dir_base}/{sample}_R1.fq.gz -2 {dir_base}/{sample}_R2.fq.gz"
    elif sample in SAMPLES_SE:
        arg = f"-U {dir_base}/{sample}.fq.gz"
    else:
        raise ValueError(f"Sample {sample} not found in SAMPLES_PE or SAMPLES_SE.")

    return arg


def _get_bowtie2_args(wildcards):
    sample = wildcards.sample
    if sample in SAMPLES_PE:
        args = "--minins 0 --maxins 2000 --no-mixed --no-discordant --dovetail"
    elif sample in SAMPLES_SE:
        args = ""
    else:
        raise ValueError(f"Sample {sample} not found in SAMPLES_PE or SAMPLES_SE.")

    return args


rule bowtie2:
    input:
        unpack(get_fastq_files),
        index=ancient(files_bowtie2_index),
    output:
        bam=protected("bowtie2/{sample}/{sample}.sorted.bam"),
        bai=protected("bowtie2/{sample}/{sample}.sorted.bam.csi"),
    log:
        "logs/{sample}/bowtie2.log",
    conda:
        "../../envs/bowtie2.yaml"
    threads: 1
    params:
        prefix=config["index_bowtie2"],
        fq=_get_fastq_input_arg,
        args=_get_bowtie2_args,
    shell:
        """
        exec 1>{log} 2>&1

        bowtie2 \\
            -p {threads} \\
            -x {params.prefix} \\
            --very-sensitive \\
            {params.args} \\
            --rg-id {wildcards.sample} \\
            --rg SM:{wildcards.sample} \\
            --rg LB:{wildcards.sample} \\
            --rg PL:ILLUMINA \\
            --rg PU:1 \\
            {params.fq} \\
                | samtools sort -@ {threads} --write-index -o {output.bam} -
        """
