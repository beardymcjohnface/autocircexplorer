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
            --novelSS \
            --readLength {params.readlen} \
            --nthread {threads} \
            --od {output.dir} \
            --tmp {params.tmpdir} \
            2> {log}
        rm -r {params.tmpdir}
        """


rule rmats_single_combine:
    input:
        rmats_out = expand(
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
        outpath = dirs["results"],
        mode = "single"
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
        g2 = ",".join(expand(os.path.join(dirs["results"], "star", "{sample}.bam"), sample=samples["group2"]))
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
        expand(
            os.path.join(dirs["results"], "rmats", "rmats_multi", "{file}"),
            file=config["rmats"]["outfiles"]["jc"] + config["rmats"]["outfiles"]["jcec"]
        )
    params:
        dir = os.path.join(dirs["results"], "rmats", "rmats_multi"),
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
            --variable-read-length \
            --novelSS \
            --nthread {threads} \
            --od {params.dir} \
            --tmp {params.tmpdir} \
            --allow-clipping \
            &> {log}
        rm -r {params.tmpdir}
        """


rule rmats_multi_combine:
    input:
        rmats_out = expand(
            os.path.join(dirs["results"], "rmats", "rmats_multi", "{file}"),
            file=config["rmats"]["outfiles"]["jc"] + config["rmats"]["outfiles"]["jcec"]),
        lib = os.path.join(dirs["results"], "lib.counts.tsv")
    output:
        long = os.path.join(dirs["results"],"rmats_multi_summary_long.csv.gz"),
        counts = expand(os.path.join(dirs["results"], "{file}.{count}.tsv"),
            file=["A3SS", "A5SS", "MXE", "RI", "SE"],
            count=["raw","CPM"])
    params:
        outpath = os.path.join(dirs["results"],"rmats", "rmats_multi"),
        respath = dirs["results"],
        mode = "multi",
        type = config["rmats"]["count"],
        s1 = samples["group1"],
        s2 = samples["group2"]
    log:
        os.path.join(dirs["logs"], "rmats_multi_combine.log")
    script:
        os.path.join(dirs["scripts"], "rmats_combine.py")


rule lib_count:
    input:
        r1=lambda wildcards: samples["reads"][wildcards.sample]["R1"],
    params:
        r2 = lambda wildcards: samples["reads"][wildcards.sample]["R2"] if samples["reads"][wildcards.sample]["R2"] else "",
        cat = lambda wildcards: "zcat" if samples["reads"][wildcards.sample]["R1"].endswith(".gz") else "cat"
    output:
        temp(os.path.join(dirs["results"], "{sample}.lib"))
    shell:
        """
        {params.cat} {input.r1} {params.r2} | wc -l | awk '{{print "{wildcards.sample}\t" $1 / 4 }}' > {output}
        """


rule combine_lib_count:
    input:
        expand(
        os.path.join(dirs["results"], "{sample}.lib"),
        sample=samples["names"]
    ),
    output:
        os.path.join(dirs["results"], "lib.counts.tsv")
    shell:
        "cat {input} > {output}"

