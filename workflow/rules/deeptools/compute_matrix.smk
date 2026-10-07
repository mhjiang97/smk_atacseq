rule generate_transcripts_bed:
    input:
        gtf=config["gtf"],
    output:
        bed="ref/transcripts.bed",
    conda:
        "../../envs/python.yaml"
    params:
        feature="transcript",
        attributes={"tag": ["Ensembl_canonical"]},
    script:
        "../../scripts/generate_transcripts_bed.py"


rule compute_matrix:
    input:
        bigwig="bowtie2/{sample}/{sample}.md.filtered.shifted.sorted.bw",
        bed="ref/transcripts.bed",
    output:
        matrix="deeptools/{sample}/{sample}.mat.gz",
    log:
        "logs/{sample}/compute_matrix.log",
    conda:
        "../../envs/deeptools.yaml"
    threads: 1
    shell:
        """
        computeMatrix \\
            scale-regions \\
            --numberOfProcessors {threads} \\
            --regionsFileName {input.bed} \\
            --scoreFileName {input.bigwig} \\
            --outFileName {output.matrix} \\
            --startLabel TSS \\
            --endLabel TES \\
            --upstream 1000 \\
            --downstream 1000 \\
            --binSize 10 \\
            --missingDataAsZero \\
            --skipZeros \\
            --samplesLabel {wildcards.sample} \\
            1>{log} 2>&1
        """


use rule compute_matrix as compute_matrix_single_end with:
    input:
        bigwig="bedtools/{sample}/{sample}.md.filtered.shifted.extended.sorted.bw",
        bed="ref/transcripts.bed",


rule compute_matrix_split:
    input:
        bigwigs=expand(
            "bowtie2/{{sample_pe}}/{{sample_pe}}.md.filtered.{type_fragment}.shifted.sorted.bw",
            type_fragment=["nfr", "1n", "2n", "3n"],
        ),
        bed="ref/tss.bed",
    output:
        matrix="deeptools/{sample_pe}/{sample_pe}.splitted.mat.gz",
    log:
        "logs/{sample_pe}/compute_matrix_split.log",
    conda:
        "../../envs/deeptools.yaml"
    threads: 1
    shell:
        """
        computeMatrix \\
            reference-point \\
            --numberOfProcessors {threads} \\
            --regionsFileName {input.bed} \\
            --scoreFileName {input.bigwigs} \\
            --outFileName {output.matrix} \\
            --referencePoint TSS \\
            --upstream 1000 \\
            --downstream 1000 \\
            --binSize 10 \\
            --missingDataAsZero \\
            --skipZeros \\
            --samplesLabel NFR Mono-N Di-N Tri-N \\
            1>{log} 2>&1
        """
