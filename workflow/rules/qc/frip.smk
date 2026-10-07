rule calculate_frip_single_end:
    input:
        bed="bedtools/{sample_se}/{sample_se}.md.filtered.nsorted.shifted.bed",
        peak="{caller}/{sample_se}/{sample_se}_peaks.narrowPeak",
    output:
        txt="frip/{sample_se}/{sample_se}.frip.{caller}.txt",
    log:
        "logs/{sample_se}/calculate_frip.{caller}.log",
    conda:
        "../../envs/bedtools.yaml"
    shell:
        """
        exec 1>{log} 2>&1

        total=$(wc -l <{input.bed})
        in_peak=$(bedtools intersect -u -a {input.bed} -b {input.peak} | wc -l)
        frip=$(echo "scale=4; ${{in_peak}} / ${{total}}" | bc)
        echo -e "Total\\tIn_peak\\tFRiP" >{output.txt}
        echo -e "${{total}}\\t${{in_peak}}\\t${{frip}}" >>{output.txt}
        """


rule calculate_frip:
    input:
        bam="bowtie2/{sample_pe}/{sample_pe}.md.filtered.shifted.sorted.bam",
        peak="{caller}/{sample_pe}/{sample_pe}_peaks.narrowPeak",
    output:
        txt="frip/{sample_pe}/{sample_pe}.frip.{caller}.txt",
    log:
        "logs/{sample_pe}/calculate_frip.{caller}.log",
    shell:
        """
        exec 1>{log} 2>&1

        total=$(samtools view -c {input.bam})
        in_peak=$(samtools view -c -L {input.peak} {input.bam})
        frip=$(echo "scale=4; ${{in_peak}} / ${{total}}" | bc)
        echo -e "Total\\tIn_peak\\tFRiP" >{output.txt}
        echo -e "${{total}}\\t${{in_peak}}\\t${{frip}}" >>{output.txt}
        """
