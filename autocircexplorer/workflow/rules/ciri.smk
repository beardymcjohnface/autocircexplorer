rule bwa_index:
    input:
        config["args"]["fa"]
    output:
        expand(
            config["args"]["fa"] + "{suffix}",
            suffix=[".amb", ".ann", ".bwt",  ".pac",  ".sa"])
    conda:
        os.path.join(dirs["envs"], "bwa.yaml")
    benchmark:
        os.path.join(dirs["bench"], "bwa_index.txt")
    log:
        os.path.join(dirs["logs"], "bwa_index.err")
    shell:
        "bwa index -a bwtsw {input} 2> {log}"


rule bwa_align:
    input:
        ref = config["args"]["fa"],
        idx = expand(
            config["args"]["fa"] + "{suffix}",
            suffix=[".amb", ".ann", ".bwt",  ".pac",  ".sa"]),
        r1=lambda wildcards: samples["reads"][wildcards.sample]["R1"],
    output:
        sam = temp(os.path.join(dirs["results"], "bwa", "{sample}.sam"))
    group:
        "ciri2"
    threads:
        config["resources"]["big"]["cpu"]
    resources:
        mem_mb = config["resources"]["big"]["mem_mb"],
        time = config["resources"]["big"]["time"]
    params:
        r2 = lambda wildcards: samples["reads"][wildcards.sample]["R2"] if samples["reads"][wildcards.sample]["R2"] else "",
    conda:
        os.path.join(dirs["envs"], "bwa.yaml")
    benchmark:
        os.path.join(dirs["bench"], "bwa_align.{sample}.txt")
    log:
        os.path.join(dirs["logs"], "bwa_align.{sample}.err")
    shell:
        "bwa mem -T 19 -t {threads} {input.ref} {input.r1} {params.r2} 2> {log} > {output}"


rule run_ciri2:
    input:
        sam = os.path.join(dirs["results"], "bwa", "{sample}.sam"),
        fa = config["args"]["fa"]
    output:
        os.path.join(dirs["results"], "ciri2", "{sample}.ciri2")
    group:
        "ciri2"
    threads:
        config["resources"]["big"]["cpu"]
    resources:
        mem_mb = config["resources"]["big"]["mem_mb"],
        time = config["resources"]["big"]["time"]
    params:
        script = os.path.join(dirs["scripts"], "CIRI2.pl")
    benchmark:
        os.path.join(dirs["bench"], "ciri2.{sample}.txt")
    log:
        os.path.join(dirs["logs"], "ciri2.{sample}.err")
    shell:
        "perl {params.script} -I {input.sam} -O {output} -F {input.fa} -T {threads} &> {log}"
