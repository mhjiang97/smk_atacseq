def _get_genrich_args(wildcards):
    sample = wildcards.sample
    if sample in SAMPLES_PE:
        return ""
    elif sample in SAMPLES_SE:
        return f"-y -d {config['extend_genrich']}"
    else:
        raise ValueError(f"Sample {sample} not found in SAMPLES_PE or SAMPLES_SE.")


rule genrich:
    input:
        bam="bowtie2/{sample}/{sample}.md.filtered.nsorted.bam",
    output:
        bed="genrich/{sample}/{sample}_peaks.narrowPeak",
    log:
        "logs/{sample}/genrich.log",
    conda:
        "../../envs/genrich.yaml"
    params:
        args=_get_genrich_args,
    shell:
        """
        Genrich \\
            {params.args} \\
            -t {input.bam} \\
            -o {output.bed} \\
            -j \\
            -p 0.01 \\
            -a 200 \\
            -l 0 \\
            -g 100 \\
            1>{log} 2>&1
        """
