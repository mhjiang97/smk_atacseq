rule bedtools_shift:
    input:
        bed="bedtools/{sample_se}/{sample_se}.md.filtered.nsorted.bed",
        fai=f"{config['fasta']}.fai",
    output:
        bed="bedtools/{sample_se}/{sample_se}.md.filtered.nsorted.shifted.bed",
    log:
        "logs/{sample_se}/bedtools_shift.log",
    conda:
        "../../envs/bedtools.yaml"
    shell:
        """
        bedtools shift \\
            -i {input.bed} \\
            -g {input.fai} \\
            -p 4 \\
            -m -5 \\
            1>{output.bed} 2>{log}
        """
