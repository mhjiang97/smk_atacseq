rule extend_bed:
    input:
        bed="bedtools/{sample_se}/{sample_se}.md.filtered.nsorted.shifted.bed",
        fai=f"{config['fasta']}.fai",
    output:
        bed="bedtools/{sample_se}/{sample_se}.md.filtered.nsorted.shifted.extended.bed",
    log:
        "logs/{sample_se}/extend_bed.log",
    params:
        _extend=config["extend_bed"],
    shell:
        """
        awk -v OFS='\\t' -v extend={params._extend} '
            FNR == NR {{ len[$1] = $2; next }}
            {{
                if ($6 == "+") {{
                    start = $2
                    end = $2 + extend
                }} else if ($6 == "-") {{
                    start = $3 - extend
                    end = $3
                }} else {{
                    next
                }}
                if (start < 0) start = 0
                if ($1 in len && end > len[$1]) end = len[$1]
                if (start < end) print $1, start, end, $4, $5, $6
            }}
        ' {input.fai} {input.bed} 1>{output.bed} 2>{log}
        """


rule bedtools_sort:
    input:
        bed="bedtools/{sample_se}/{sample_se}.md.filtered.nsorted.shifted.extended.bed",
        fai=f"{config['fasta']}.fai",
    output:
        bed="bedtools/{sample_se}/{sample_se}.md.filtered.shifted.extended.sorted.bed",
    log:
        "logs/{sample_se}/bedtools_sort.log",
    conda:
        "../../envs/bedtools.yaml"
    shell:
        """
        bedtools sort \\
            -faidx {input.fai} \\
            -i {input.bed} \\
            1>{output.bed} 2>{log}
        """
