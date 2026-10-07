rule macs3_paired_end:
    input:
        bam="bowtie2/{sample_pe}/{sample_pe}.md.filtered.shifted.sorted.bam",
    output:
        bed="macs3/{sample_pe}/{sample_pe}_peaks.narrowPeak",
    log:
        "logs/{sample_pe}/macs3_paired_end.log",
    conda:
        "../../envs/macs3.yaml"
    params:
        prefix="macs3/{sample_pe}/{sample_pe}",
        size_genome=SIZES_GENOME[config["genome"]],
        q=0.01,
    shell:
        """
        macs3 callpeak \\
            --format BAMPE \\
            --bdg \\
            --gsize {params.size_genome} \\
            --qvalue {params.q} \\
            --keep-dup all \\
            --treatment {input.bam} \\
            --name {params.prefix} \\
            1>{log} 2>&1
        """


rule macs3_single_end:
    input:
        bed="bedtools/{sample_se}/{sample_se}.md.filtered.nsorted.shifted.bed",
    output:
        bed="macs3/{sample_se}/{sample_se}_peaks.narrowPeak",
    log:
        "logs/{sample_se}/macs3_single_end.log",
    conda:
        "../../envs/macs3.yaml"
    params:
        prefix="macs3/{sample_se}/{sample_se}",
        size_genome=SIZES_GENOME[config["genome"]],
        q=0.01,
        extsize=config["extend_macs3"],
        shift=int(-1 * config["extend_macs3"] / 2),
    shell:
        """
        macs3 callpeak \\
            --format BED \\
            --bdg \\
            --nomodel \\
            --shift {params.shift} \\
            --extsize {params.extsize} \\
            --gsize {params.size_genome} \\
            --qvalue {params.q} \\
            --keep-dup all \\
            --treatment {input.bed} \\
            --name {params.prefix} \\
            1>{log} 2>&1
        """
