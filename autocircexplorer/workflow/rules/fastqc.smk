import os.path

rule fastqc:
    input:
        lambda wildcards: samples["reads"][wildcards.sample][wildcards.read]
    output:
        html1=os.path.join(dirs["results"], "multiqc", "{sample}_{read}_fastqc.html"),
        zip1=os.path.join(dirs["results"], "multiqc", "{sample}_{read}_fastqc.zip")
    params:
        extra = "--quiet",
        tmp = os.path.join(dirs["results"], "multiqc", "{sample}_{read}.fastq.gz"),
        dir = os.path.join(dirs["results"], "multiqc")
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


rule run_multiqc:
    input:
        targets["fastqc"],
        targets["flagstat"]
    output:
        html = os.path.join(dirs["results"], "multiqc_report.html"),
        dir = directory(os.path.join(dirs["results"], "multiqc_report_data"))
    params:
        os.path.join(dirs["results"], "multiqc")
    log:
        os.path.join(dirs["logs"], "run_multiqc.log")
    benchmark:
        os.path.join(dirs["bench"], "run_multiqc.txt")
    conda:
        os.path.join(dirs["envs"], "multiqc.yaml")
    shell:
        "multiqc {params} --filename {output.html} 2> {log}"
