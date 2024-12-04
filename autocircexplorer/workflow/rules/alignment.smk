rule star_align:
    input:
        ref = config["args"]["ref"],
        r1=lambda wildcards: samples["reads"][wildcards.sample]["R1"],
    output:
        bam = os.path.join(dirs["results"], "star", "{sample}.bam"),
        out = expand(
            os.path.join(dirs["results"], "star", "{{sample}}.{file}"),
            file=[
                "Chimeric.out.junction",
                "Log.out",
                "Log.final.out",
                "Log.progress.out",
                "Log.std.out",
                "SJ.out.tab"]
            )
    threads:
        config["resources"]["big"]["cpu"]
    resources:
        mem_mb = config["resources"]["big"]["mem_mb"],
        time = config["resources"]["big"]["time"]
    params:
        star = config["starparams"],
        prefix = os.path.join(dirs["results"], "star", "{sample}."),
        r2 = lambda wildcards: samples["reads"][wildcards.sample]["R2"] if samples["reads"][wildcards.sample]["R2"] else "",
        readcmd = lambda wildcards: "zcat" if samples["reads"][wildcards.sample]["R1"].endswith(".gz") else "cat"
    conda:
        os.path.join(dirs["envs"], "star.yaml")
    benchmark:
        os.path.join(dirs["bench"], "star_align.{sample}.txt")
    log:
        os.path.join(dirs["logs"], "star_align.{sample}.err")
    shell:
        """
        STAR {params.star} \
            --runThreadN {threads} \
            --genomeDir {input.ref} \
            --readFilesCommand {params.readcmd} \
            --outStd SAM \
            --outFileNamePrefix {params.prefix} \
            --readFilesIn {input.r1} {params.r2} \
        | samtools sort -@ {threads} -m 1G \
        > {output.bam}
        """


rule bam_index:
    input:
        "{filepath}.bam"
    output:
        "{filepath}.bam.bai"
    conda:
        os.path.join(dirs["envs"], "star.yaml")
    shell:
        "samtools index {input} -o {output}"

