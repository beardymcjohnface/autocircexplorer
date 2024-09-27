rule rmats_single_sample_file:
    input:
        os.path.join(dirs["results"], "star", "{sample}.bam")
    output:
        os.path.join(dirs["results"], "rmats", "{sample}.b1")
    shell:
        """echo {input} > {output}"""


rule rmats_single:
    input:
        bam = os.path.join(dirs["results"], "star", "{sample}.bam"),
        b1 = os.path.join(dirs["results"], "rmats", "{sample}.b1"),
        gtf = config["args"]["gtf"]
    output:
        files = expand(
            os.path.join(dirs["results"], "rmats", "{{sample}}.rmats_single", "{file}"),
            file=config["rmats"]["outfiles"]["jc"] + config["rmats"]["outfiles"]["jcec"]
        ),
        dir = directory(os.path.join(dirs["results"], "rmats", "{sample}.rmats_single"))
    params:
        paired = lambda w: "paired" if samples["reads"][w.sample]["R2"] else "single",
        readlen = config["args"]["readlen"],
        tmpdir = os.path.join(dirs["results"], "rmats", "{sample}.rmats_single_tmp")
    threads:
        8
    singularity:
        "docker://xinglab/rmats"
    log:
        os.path.join(dirs["logs"], "rmats_single.{sample}.log")
    shell:
        """
        python /rmats/rmats.py \
            --b1 {input.b1} \
            --statoff \
            --gtf {input.gtf} \
            -t {params.paired} \
            --readLength {params.readlen} \
            --nthread {threads} \
            --od {output.dir} \
            --tmp {params.tmpdir} \
            2> {log}
        rm -r {params.tmpdir}
        """


rule rmats_single_combine:
    input:
        expand(
            os.path.join(dirs["results"], "rmats", "{sample}.rmats_single", "{file}"),
            sample=samples["names"],
            file=config["rmats"]["outfiles"]["jc"] + config["rmats"]["outfiles"]["jcec"]
        )
    output:
        long = os.path.join(dirs["results"], "rmats_single_summary_long.csv.gz"),
        counts = expand(
            os.path.join(dirs["results"], "{event}.{count}.csv.gz"),
            event=["A3SS","A5SS","MXE","RI","SE"],
            count=["skip","inclusion"]
        )
    params:
        outpath = dirs["results"]
    log:
        os.path.join(dirs["logs"], "rmats_single_combine.log")
    script:
        os.path.join(dirs["scripts"], "rmats_combine.py")



rule rmats_multi_sample_files:
    input:
        expand(os.path.join(dirs["results"], "star", "{sample}.bam"), sample=samples["names"])
    output:
        g1 = os.path.join(dirs["results"], "rmats", "group1.txt"),
        g2 = os.path.join(dirs["results"], "rmats", "group2.txt")
    params:
        g1 = ",".join(expand(os.path.join(dirs["results"], "star", "{sample}.bam"), sample=samples["group1"])),
        g2 = ",".join(expand(os.path.join(dirs["results"], "star", "{sample}.bam"), sample=samples["group2"])),
    shell:
        """
        echo {params.g1} > {output.g1}
        echo {params.g2} > {output.g2}
        """


rule rmats_multi:
    input:
        bams = expand(os.path.join(dirs["results"], "star", "{sample}.bam"), sample=samples["names"]),
        b1 = os.path.join(dirs["results"], "rmats", "group1.txt"),
        b2 = os.path.join(dirs["results"], "rmats", "group2.txt"),
        gtf = config["args"]["gtf"]
    output:
        directory(os.path.join(dirs["results"], "rmats", "rmats_multi"))
    params:
        paired = lambda w: "paired" if samples["reads"][samples["names"][0]]["R2"] else "single",
        readlen = config["args"]["readlen"],
        tmpdir = os.path.join(dirs["results"], "rmats", "rmats_multi_tmp")
    threads:
        40
    singularity:
        "docker://xinglab/rmats"
    log:
        os.path.join(dirs["logs"], "rmats_multi.log")
    shell:
        """
        python /rmats/rmats.py \
            --b1 {input.b1} \
            --b2 {input.b2} \
            --gtf {input.gtf} \
            -t {params.paired} \
            --readLength {params.readlen} \
            --nthread {threads} \
            --od {output} \
            --tmp {params.tmpdir} \
            2> {log}
        rm -r {params.tmpdir}
        """
