rule samtools_sort:
    input:
        bam="bowtie2/{sample}/{sample}.sorted.md.filtered.bam",
    output:
        bam="bowtie2/{sample}/{sample}.md.filtered.nsorted.bam",
    log:
        "logs/{sample}/samtools_sort.log",
    threads: 1
    shell:
        """
        samtools sort -@ {threads} -n -o {output.bam} {input.bam} 1>{log} 2>&1
        """
