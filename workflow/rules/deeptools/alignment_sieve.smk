rule alignment_sieve:
    input:
        bam="bowtie2/{sample_pe}/{sample_pe}.sorted.md.filtered.bam",
    output:
        tmp=temp("bowtie2/{sample_pe}/shifted.tmp.bam"),
        bam="bowtie2/{sample_pe}/{sample_pe}.md.filtered.shifted.sorted.bam",
    log:
        "logs/{sample_pe}/alignment_sieve.log",
    conda:
        "../../envs/deeptools3.yaml"
    threads: 1
    shell:
        """
        exec 1>{log} 2>&1

        alignmentSieve \\
            --numberOfProcessors {threads} \\
            --bam {input.bam} \\
            --outFile {output.tmp} \\
            --ATACshift

        samtools sort -@ {threads} --write-index -o {output.bam} {output.tmp}
        """


rule alignment_sieve_split:
    input:
        bam="bowtie2/{sample_pe}/{sample_pe}.sorted.md.filtered.{type_fragment}.bam",
    output:
        bam="bowtie2/{sample_pe}/{sample_pe}.md.filtered.{type_fragment}.shifted.sorted.bam",
        tmp=temp("bowtie2/{sample_pe}/shifted.tmp.{type_fragment}.bam"),
    log:
        "logs/{sample_pe}/alignment_sieve_split.{type_fragment}.log",
    conda:
        "../../envs/deeptools3.yaml"
    threads: 1
    shell:
        """
        exec 1>{log} 2>&1

        alignmentSieve \\
            --numberOfProcessors {threads} \\
            --bam {input.bam} \\
            --outFile {output.tmp} \\
            --ATACshift

        samtools sort -@ {threads} --write-index -o {output.bam} {output.tmp}
        """
