_pdfs = multiext(
    "bowtie2/{sample}/{sample}.sorted.md.filtered",
    ".base_distribution_by_cycle.pdf",
    ".quality_by_cycle.pdf",
    ".quality_distribution.pdf",
    ".read_length_histogram.pdf",
)


rule collect_multiple_metrics:
    input:
        bam="bowtie2/{sample}/{sample}.sorted.md.filtered.bam",
        fasta=config["fasta"],
    output:
        pdfs=_pdfs,
    log:
        "logs/{sample}/collect_multiple_metrics.log",
    conda:
        "../../envs/gatk.yaml"
    resources:
        mem_mb=1,
        tmpdir=lambda wildcards: f"bowtie2/{wildcards.sample}",
    params:
        prefix="bowtie2/{sample}/{sample}.sorted.md.filtered",
    shell:
        """
        gatk CollectMultipleMetrics \\
            --java-options \"-Xmx{resources.mem_mb}M -XX:-UsePerfData\" \\
            --INPUT {input.bam} \\
            --OUTPUT {params.prefix} \\
            --REFERENCE_SEQUENCE {input.fasta} \\
            --TMP_DIR {resources.tmpdir} \\
            1>{log} 2>&1
        """
