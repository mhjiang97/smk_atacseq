def _cpm_scale(wildcards, input):
    with open(input.bed) as fh:
        n_reads = sum(1 for line in fh if line.strip())

    return 1_000_000 / n_reads


rule bedtools_genomecov:
    input:
        bed="bedtools/{sample_se}/{sample_se}.md.filtered.shifted.extended.sorted.bed",
        fai=f"{config['fasta']}.fai",
    output:
        bedgraph="bedtools/{sample_se}/{sample_se}.md.filtered.shifted.extended.sorted.bedGraph",
    log:
        "logs/{sample_se}/bedtools_genomecov.log",
    conda:
        "../../envs/bedtools.yaml"
    params:
        scale=_cpm_scale,
    shell:
        """
        bedtools genomecov \\
            -bg \\
            -scale {params.scale} \\
            -i {input.bed} \\
            -g {input.fai} \\
            1>{output.bedgraph} 2>{log}
        """
