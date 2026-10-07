rule plot_heatmap:
    input:
        matrix="deeptools/{sample}/{sample}.mat.gz",
    output:
        pdf="deeptools/{sample}/{sample}.heatmap.pdf",
    log:
        "logs/{sample}/plot_heatmap.log",
    conda:
        "../../envs/deeptools.yaml"
    shell:
        """
        plotHeatmap \\
            --matrixFile {input.matrix} \\
            --outFileName {output.pdf} \\
            1>{log} 2>&1
        """


rule plot_heatmap_split:
    input:
        matrix="deeptools/{sample_pe}/{sample_pe}.splitted.mat.gz",
    output:
        pdf="deeptools/{sample_pe}/{sample_pe}.splitted.heatmap.pdf",
    log:
        "logs/{sample_pe}/plot_heatmap_split.log",
    conda:
        "../../envs/deeptools.yaml"
    shell:
        """
        plotHeatmap \\
            --matrixFile {input.matrix} \\
            --outFileName {output.pdf} \\
            1>{log} 2>&1
        """
