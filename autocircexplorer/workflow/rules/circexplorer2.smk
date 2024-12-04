rule parse_star:
    input:
        os.path.join(dirs["results"], "star", "{sample}.Chimeric.out.junction")
    output:
        os.path.join(dirs["results"], "ce2", "{sample}.circexplorer2.parse")
    conda:
        os.path.join(dirs["envs"], "ce2.yaml")
    benchmark:
        os.path.join(dirs["bench"], "parse_star.{sample}.txt")
    log:
        os.path.join(dirs["logs"], "parse_star.{sample}.err")
    shell:
        "CIRCexplorer2 parse -t STAR {input} -b {output} > {log}"


rule gtf_to_annotation:
    input:
        config["args"]["gtf"]
    output:
        os.path.join(dirs["results"], "ce2", "refFlatFile.txt")
    benchmark:
        os.path.join(dirs["bench"], "gtf_to_annotation.txt")
    log:
        os.path.join(dirs["logs"], "gtf_to_annotation.err")
    script:
        os.path.join(dirs["scripts"], "gtf_to_flatfile.py")


rule circexplorer2_annotate:
    input:
        flat = os.path.join(dirs["results"], "ce2", "refFlatFile.txt"),
        fa = config["args"]["fa"],
        parse = os.path.join(dirs["results"], "ce2", "{sample}.circexplorer2.parse")
    output:
        os.path.join(dirs["results"], "ce2", "{sample}.circexplorer2.annotated")
    conda:
        os.path.join(dirs["envs"], "ce2.yaml")
    benchmark:
        os.path.join(dirs["bench"], "circexplorer2_annotate.{sample}.txt")
    log:
        os.path.join(dirs["logs"], "circexplorer2_annotate.{sample}.err")
    shell:
        "CIRCexplorer2 annotate "
            "-r {input.flat} "
            "-g {input.fa} "
            "-b {input.parse} "
            "-o {output} "
            "&> {log}; "
