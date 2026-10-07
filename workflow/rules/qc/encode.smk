rule tss_enrich:
    input:
        bam="bowtie2/{sample}/{sample}.sorted.md.filtered.bam",
        tss="ref/tss.bed",
        fai=f"{config['fasta']}.fai",
    output:
        pdf="encode/{sample}/{sample}.tss_enrich.pdf",
        pdf_large="encode/{sample}/{sample}.large_tss_enrich.pdf",
        qc="encode/{sample}/{sample}.tss_enrich.qc",
    log:
        "logs/{sample}/tss_enrich.log",
    conda:
        "../../envs/python.yaml"
    threads: 1
    params:
        len_read=None,
    script:
        "../../scripts/tss_enrich.py"
