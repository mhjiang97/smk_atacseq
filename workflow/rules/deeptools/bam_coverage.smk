rule bam_coverage:
    input:
        bam="bowtie2/{sample_pe}/{sample_pe}.md.filtered.shifted.sorted.bam",
    output:
        bigwig="bowtie2/{sample_pe}/{sample_pe}.md.filtered.shifted.sorted.bw",
    log:
        "logs/{sample_pe}/bam_coverage.log",
    conda:
        "../../envs/deeptools.yaml"
    threads: 1
    params:
        mode_normalization=config["mode_normalization"],
        size_bin=1,
        size_genome=SIZES_GENOME[config["genome"]],
    shell:
        """
        bamCoverage \\
            --numberOfProcessors {threads} \\
            --bam {input.bam} \\
            --outFileFormat bigwig \\
            --outFileName {output.bigwig} \\
            --binSize {params.size_bin} \\
            --normalizeUsing {params.mode_normalization} \\
            --effectiveGenomeSize {params.size_genome} \\
            --skipNAs \\
            --extendReads \\
            1>{log} 2>&1
        """


rule bam_coverage_split:
    input:
        bam="bowtie2/{sample_pe}/{sample_pe}.md.filtered.{type_fragment}.shifted.sorted.bam",
    output:
        bigwig="bowtie2/{sample_pe}/{sample_pe}.md.filtered.{type_fragment}.shifted.sorted.bw",
    log:
        "logs/{sample_pe}/bam_coverage_split.{type_fragment}.log",
    conda:
        "../../envs/deeptools.yaml"
    threads: 1
    params:
        mode_normalization=config["mode_normalization"],
        size_bin=1,
        size_genome=SIZES_GENOME[config["genome"]],
    shell:
        """
        bamCoverage \\
            --numberOfProcessors {threads} \\
            --bam {input.bam} \\
            --outFileFormat bigwig \\
            --outFileName {output.bigwig} \\
            --binSize {params.size_bin} \\
            --normalizeUsing {params.mode_normalization} \\
            --effectiveGenomeSize {params.size_genome} \\
            --skipNAs \\
            --extendReads \\
            1>{log} 2>&1
        """
