rule samtools_faidx:
    input:
        fasta=config["fasta"],
    output:
        fai=f"{config['fasta']}.fai",
    log:
        "logs/samtools_faidx.log",
    shell:
        """
        samtools faidx {input.fasta} --fai-idx {output.fai} 1>{log} 2>&1
        """
