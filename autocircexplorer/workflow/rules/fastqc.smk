import os.path

rule fastqc:
    input:
        lambda wildcards: samples["reads"][wildcards.sample][wildcards.read]
    output:
        html1=os.path.join(dirs["results"], "fastqc", "{sample}_{read}_fastqc.html"),
        zip1=os.path.join(dirs["results"], "fastqc", "{sample}_{read}_fastqc.zip")
    params:
        extra = "--quiet",
        tmp = os.path.join(dirs["results"], "fastqc", "{sample}_{read}.fastq.gz"),
        dir = os.path.join(dirs["results"], "fastqc")
    log:
        os.path.join(dirs["logs"], "fastqc.{sample}_{read}.log")
    benchmark:
        os.path.join(dirs["bench"], "fastqc.{sample}_{read}.txt")
    conda:
        os.path.join(dirs["envs"], "fastqc.yaml")
    threads: 1
    resources:
        mem_mb = 1024
    shell:
        "ln -sr {input} {params.tmp}; "
        "fastqc {params.tmp} -t {threads} --outdir {params.dir} 2> {log}; "


rule multiqc_fastqc:
    input:
        targets["fastqc"]
    output:
        os.path.join(dirs["results"], "multi_fastqc_report.html")
    params:
        dir = os.path.join(dirs["results"], "fastqc")
    log:
        os.path.join(dirs["logs"], "multiqc_fastqc.log")
    benchmark:
        os.path.join(dirs["bench"], "multiqc_fastqc.txt")
    conda:
        os.path.join(dirs["envs"], "multiqc.yaml")
    shell:
        "multiqc {params.dir} --filename {output} 2> {log}"
