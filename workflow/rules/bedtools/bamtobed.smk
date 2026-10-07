rule bedtools_bamtobed:
    input:
        bam="bowtie2/{sample_se}/{sample_se}.md.filtered.nsorted.bam",
    output:
        bed="bedtools/{sample_se}/{sample_se}.md.filtered.nsorted.bed",
    log:
        "logs/{sample_se}/bedtools_bamtobed.log",
    conda:
        "../../envs/bedtools.yaml"
    shell:
        """
        bedtools bamtobed -i {input.bam} 1>{output.bed} 2>{log}
        """
