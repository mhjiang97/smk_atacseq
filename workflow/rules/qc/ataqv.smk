rule generate_tss_bed:
    input:
        gtf=config["gtf"],
    output:
        bed="ref/tss.bed",
    params:
        feature="transcript",
        attributes={"tag": ["Ensembl_canonical"]},
    run:
        idx_chrom = 0
        idx_feature = 2
        idx_start = 3
        idx_end = 4
        idx_strand = 6
        idx_attributes = 8

        def parse_attributes(text):
            parsed = {}
            for item in text.strip().rstrip(";").split(";"):
                item = item.strip()
                if not item:
                    continue
                key, _, rest = item.partition(" ")
                parsed.setdefault(key, set()).add(rest.strip().strip('"'))
            return parsed

        def matches_attributes(text, required):
            parsed = parse_attributes(text)
            return all(
                parsed.get(key, set()) & set(values)
                for key, values in required.items()
            )

        tss = set()
        with open(input.gtf, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if fields[idx_feature] != params.feature:
                    continue
                if params.attributes and not matches_attributes(
                    fields[idx_attributes], params.attributes
                ):
                    continue
                start = int(fields[idx_start])
                end = int(fields[idx_end])
                strand = fields[idx_strand]
                if strand == "+":
                    start_bed = start - 1
                    end_bed = start
                else:
                    start_bed = end - 1
                    end_bed = end
                chrom = fields[idx_chrom]
                tss.add((chrom, start_bed, end_bed, strand))
        with open(output.bed, "w") as f:
            for chrom, start, end, strand in sorted(tss):
                f.write(f"{chrom}\t{start}\t{end}\t.\t0\t{strand}\n")


rule ataqv:
    input:
        bam="bowtie2/{sample}/{sample}.sorted.md.filtered.bam",
        peak="{caller}/{sample}/{sample}_peaks.narrowPeak",
        tss="ref/tss.bed",
    output:
        metrics="ataqv/{sample}/ataqv.{caller}.json",
        ataqv="ataqv/{sample}/ataqv.{caller}.out",
    log:
        "logs/{sample}/ataqv.{caller}.log",
    conda:
        "../../envs/ataqv.yaml"
    threads: 1
    params:
        organism=config["organism"],
    shell:
        """
        ataqv \\
            --threads {threads} \\
            --peak-file {input.peak} \\
            --name {wildcards.sample} \\
            --metrics-file {output.metrics} \\
            --tss-file {input.tss} \\
            --ignore-read-groups \\
            {params.organism} \\
            {input.bam} \\
            1>{output.ataqv} 2>{log}
        """


rule mkarv:
    input:
        expand("ataqv/{sample}/ataqv.{{caller}}.json", sample=SAMPLES_PE),
    output:
        _dir=directory("ataqv/mkarv.{caller}"),
        html="ataqv/mkarv.{caller}/index.html",
    log:
        "logs/mkarv.{caller}.log",
    conda:
        "../../envs/ataqv.yaml"
    threads: 1
    shell:
        """
        mkarv --concurrency {threads} --force {output._dir} {input} 1>{log} 2>&1
        """
