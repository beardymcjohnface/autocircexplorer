import os.path

rule salmon_alignment:
    input:
        ref = config["args"]["ref"],
        r1 = lambda wildcards: samples["reads"][wildcards.sample]["R1"],
    output:
        os.path.join(dirs["results"],"salmon","{sample}","quant.sf")
    params:
        dir = os.path.join(dirs["results"],"salmon","{sample}"),
        r2 = lambda wildcards: "-2 " + samples["reads"][wildcards.sample]["R2"] if samples["reads"][wildcards.sample]["R2"] else "",
        r1 = lambda wildcards: "-1 " if samples["reads"][wildcards.sample]["R2"] else "-r "
    threads:
        config["resources"]["big"]["cpu"]
    resources:
        mem_mb = config["resources"]["big"]["mem_mb"],
        time = config["resources"]["big"]["time"]
    conda:
        os.path.join(dirs["envs"], "salmon.yaml")
    benchmark:
        os.path.join(dirs["bench"], "salmon_alignment.{sample}.txt")
    log:
        os.path.join(dirs["logs"], "salmon_alignment.{sample}.err")
    shell:
        "salmon quant "
            "-i {input.ref} "
            "-l A "
            "{params.r1} {input.r1} "
            "{params.r2} "
            "-p {threads} "
            "-o {params.dir} "
            "2> {log}"


rule combine_salmon:
    input:
        expand(os.path.join(dirs["results"],"salmon","{sample}","quant.sf"), sample=samples["names"])
    output:
        tpm = os.path.join(dirs["results"],"salmon","TPM.tsv"),
        num = os.path.join(dirs["results"],"salmon","NumReads.tsv")
    benchmark:
        os.path.join(dirs["bench"], "combine_salmon.txt")
    params:
        samples = samples["names"],
        dir = os.path.join(dirs["results"],"salmon"),
        file = "quant.sf"
    log:
        os.path.join(dirs["logs"], "combine_salmon.err")
    script:
        os.path.join(dirs["scripts"], "combine_salmon.py")


