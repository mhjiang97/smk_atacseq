def _get_fastp_io_args(wildcards, input, output):
    sample = wildcards.sample
    if sample in SAMPLES_PE or sample in CONTROLS_PE:
        return (
            f"--detect_adapter_for_pe "
            f"--in1 {input.fq_1} --in2 {input.fq_2} "
            f"--out1 {output.fq_1} --out2 {output.fq_2} "
            f"--unpaired1 {output.fq_1_unpaired} --unpaired2 {output.fq_2_unpaired}"
        )

    return f"-i {input.fq} -o {output.fq}"


rule fastp:
    input:
        fq1=f"{DIR_DATA}/{{sample}}_R1.fq.gz",
        fq2=f"{DIR_DATA}/{{sample}}_R2.fq.gz",
    output:
        fq1="fastp/{sample}/{sample}_R1.fq.gz",
        fq2="fastp/{sample}/{sample}_R2.fq.gz",
        fq1_unpaired="fastp/{sample}/{sample}.unpaired_R1.fq.gz",
        fq2_unpaired="fastp/{sample}/{sample}.unpaired_R2.fq.gz",
        fq_failed="fastp/{sample}/{sample}.failed.fq.gz",
        json="fastp/{sample}/{sample}.json",
        html="fastp/{sample}/{sample}.html",
    log:
        "logs/{sample}/fastp.log",
    threads: 1
    params:
        args=get_extra_arguments("fastp"),
        io=_get_fastp_io_args,
    shell:
        """
        fastp \\
            {params.args} \\
            --thread {threads} \\
            {params.io} \\
            --failed_out {output.fq_failed} \\
            --json {output.json} \\
            --html {output.html} \\
            1>{log} 2>&1
        """


use rule fastp as fastp_single_end with:
    input:
        fq=f"{DIR_DATA}/{{sample}}.fq.gz",
    output:
        fq=f"fastp/{{sample}}/{{sample}}.fq.gz",
        fq_failed="fastp/{sample}/{sample}.failed.fq.gz",
        json="fastp/{sample}/{sample}.json",
        html="fastp/{sample}/{sample}.html",
