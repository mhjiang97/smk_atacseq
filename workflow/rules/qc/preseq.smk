rule preseq_lc_extrap:
    input:
        bam="bowtie2/{sample}/{sample}.sorted.md.filtered.bam",
    output:
        txt="bowtie2/{sample}/{sample}.sorted.md.filtered.lc_extrap.txt",
    log:
        "logs/{sample}/preseq.lc_extrap.log",
    shadow:
        "minimal"
    conda:
        "../../envs/preseq.yaml"
    shell:
        """
        preseq lc_extrap -v -B -o {output.txt} {input.bam} 1>{log} 2>&1
        """


rule preseq_c_curve:
    input:
        bam="bowtie2/{sample}/{sample}.sorted.md.filtered.bam",
    output:
        txt="bowtie2/{sample}/{sample}.sorted.md.filtered.c_curve.txt",
    log:
        "logs/{sample}/preseq.c_curve.log",
    shadow:
        "minimal"
    conda:
        "../../envs/preseq.yaml"
    shell:
        """
        preseq c_curve -v -B -P -o {output.txt} {input.bam} 1>{log} 2>&1
        """


rule preseq_bound_pop:
    input:
        bam="bowtie2/{sample}/{sample}.sorted.md.filtered.bam",
    output:
        txt="bowtie2/{sample}/{sample}.sorted.md.filtered.bound_pop.txt",
    log:
        "logs/{sample}/preseq.bound_pop.log",
    shadow:
        "minimal"
    conda:
        "../../envs/preseq.yaml"
    shell:
        """
        preseq bound_pop -v -B -P -o {output.txt} {input.bam} 1>{log} 2>&1
        """
