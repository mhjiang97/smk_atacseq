rule bedgraph_to_bigwig:
    input:
        bedgraph="bedtools/{sample_se}/{sample_se}.md.filtered.shifted.extended.sorted.bedGraph",
        fai=f"{config['fasta']}.fai",
    output:
        bigwig="bedtools/{sample_se}/{sample_se}.md.filtered.shifted.extended.sorted.bw",
    log:
        "logs/{sample_se}/bedgraph_to_bigwig.log",
    conda:
        "../../envs/ucsc.yaml"
    shell:
        """
        bedGraphToBigWig \\
            {input.bedgraph} \\
            {input.fai} \\
            {output.bigwig} \\
            1>{log} 2>&1
        """
