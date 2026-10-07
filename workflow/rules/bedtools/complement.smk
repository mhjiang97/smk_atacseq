rule bedtools_sort_blacklist:
    input:
        blacklist=config["blacklist"],
        fai=f"{config['fasta']}.fai",
    output:
        bed="ref/blacklist.sorted.bed",
    log:
        "logs/bedtools_sort_blacklist.log",
    conda:
        "../../envs/bedtools.yaml"
    shell:
        """
        bedtools sort \\
            -faidx {input.fai} \\
            -i {input.blacklist} \\
            1>{output.bed} 2>{log}
        """


rule bedtools_complement:
    input:
        bed="ref/blacklist.sorted.bed",
        fai=f"{config['fasta']}.fai",
    output:
        bed="ref/include.bed",
    log:
        "logs/bedtools_complement.log",
    conda:
        "../../envs/bedtools.yaml"
    shell:
        """
        bedtools complement \\
            -i {input.bed} \\
            -g {input.fai} \\
            1>{output.bed} 2>{log}
        """
