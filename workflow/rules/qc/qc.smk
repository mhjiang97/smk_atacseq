rule fastqc:
    input:
        lambda wildcards: get_fastq_files(wildcards, named=False),
    output:
        _dir=branch(
            TO_CLEAN_FQ,
            then=directory("fastqc/fastp/{sample}"),
            otherwise=directory("fastqc/{sample}"),
        ),
    log:
        "logs/{sample}/fastqc.log",
    threads: 1
    shell:
        """
        mkdir -p {output._dir}
        fastqc -t {threads} -o {output._dir} {input} 1>{log} 2>&1
        """


rule multiqc:
    input:
        fastqcs=branch(
            TO_CLEAN_FQ,
            then=expand("fastqc/fastp/{sample}", sample=SAMPLES),
            otherwise=expand("fastqc/{sample}", sample=SAMPLES),
        ),
    output:
        _dir=branch(
            TO_CLEAN_FQ,
            then=directory("multiqc/fastp"),
            otherwise=directory("multiqc"),
        ),
        html=branch(
            TO_CLEAN_FQ,
            then="multiqc/fastp/multiqc_report.html",
            otherwise="multiqc/multiqc_report.html",
        ),
    log:
        "logs/multiqc.log",
    shell:
        """
        multiqc -o {output._dir} --force {input.fastqcs} 1>{log} 2>&1
        """
