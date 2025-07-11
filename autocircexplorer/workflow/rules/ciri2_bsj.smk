# pull trimmed read IDs for each BSJ and reassemble
    # build junction target list from result file?

rule ciri2_sample_bsj_reads:
    input:
        ciri2 = os.path.join(dirs["results"],"ciri2","{sample}.ciri2"),
        r1= lambda wildcards: samples["reads"][wildcards.sample]["R1"]
    output:
        os.path.join(dirs["results"],"ciri2","{sample}.reads.pkl")
    params:
        r2= lambda wildcards: samples["reads"][wildcards.sample]["R2"] if samples["reads"][wildcards.sample]["R2"] else "",
        min_reads = config["ciri2"]["min_bsj_reads"]
    log:
        os.path.join(dirs["logs"], "ciri2_sample_bsj_reads.{sample}.log")
    script:
        os.path.join(dirs["scripts"],"ciri2_sample_bsj_reads.py")


rule ciri2_combine_bsj_reads:
    input:
        expand(os.path.join(dirs["results"],"ciri2","{sample}.reads.pkl"), sample=samples["names"])
    output:
        os.path.join(dirs["results"], "ciri2", "bsj_reads.R1.fastq")
    params:
        os.path.join(dirs["results"], "ciri2", "bsj_reads.R2.fastq")
    run:
        import pickle
        r1_fh = open(output[0], "w")
        r2_fh = open(params[0], "w")
        for pkl_file in input:
            with open(pkl_file, "rb") as fh:
                sample_reads = pickle.load(fh)
                for bsj_id in sample_reads["R1"].keys():
                    for line in sample_reads["R1"][bsj_id]:
                        r1_fh.write(line)
                    if bsj_id in sample_reads["R2"]:
                        for line in sample_reads["R2"][bsj_id]:
                            r2_fh.write(line)
        r1_fh.close()
        r2_fh.close()


rule ciri2_fastp_bsj_reads:
    input:
        os.path.join(dirs["results"],"ciri2","bsj_reads.R1.fastq")
    output:
        os.path.join(dirs["results"],"ciri2","bsj_reads.R1.trimmed.fastq")
    params:
        inr2 = os.path.join(dirs["results"], "ciri2", "bsj_reads.R2.fastq"),
        outr2 = os.path.join(dirs["results"], "ciri2", "bsj_reads.R2.trimmed.fastq")
    conda:
        os.path.join(dirs["envs"], "fastp.yaml")
    shell:
        """
        if [ -s {params.inr2} ]
            then fastp -i {input} -I {params.inr2} -o {output} -O {params.outr2}
            else fastp -i {input} -o {output}
        fi
        """


rule assemble_ciri2_bsj:
    input:
        os.path.join(dirs["results"],"ciri2","bsj_reads.R1.trimmed.fastq")
    output:
        os.path.join(dirs["results"],"ciri2","bsj_assembly.fasta")
    params:
        r2 = os.path.join(dirs["results"],"ciri2","bsj_reads.R2.trimmed.fastq"),
        dir = os.path.join(dirs["results"],"ciri2","trinity_out"),
        tmp = os.path.join(dirs["results"],"ciri2","trinity_out.Trinity.fasta")
    threads:
        20
    resources:
        mem = "128G"
    conda:
        os.path.join(dirs["envs"], "trinity.yaml")
    shell:
        """
        if [ -s {params.r2} ]
            then Trinity \
                --seqType fq \
                --max_memory {resources.mem} \
                --left {input}  \
                --right {params.r2} \
                --CPU {threads} \
                --output {params.dir} \
                --full_cleanup
            else Trinity \
                --seqType fq \
                --max_memory {resources.mem} \
                --single {input}  \
                --CPU {threads} \
                --output {params.dir} \
                --full_cleanup
        fi
        mv {params.tmp} {output}
        """


rule ciri2_reference_bsj_seqs:
    input:
        ciri2 = expand(os.path.join(dirs["results"],"ciri2","{sample}.ciri2"), sample=samples["names"]),
        ref = config["args"]["fa"]
    output:
        os.path.join(dirs["results"], "ciri2", "ref_bsj_seqs.fa")
    conda:
        os.path.join(dirs["envs"], "star.yaml")
    shell:
        """
        cat {input.ciri2} \
            | cut -f2,3,4 \
            | sort \
            | uniq \
            | awk '{{print ">"$1":"$2"|"$3; system("samtools faidx {input.ref} " $1":"$2"-"$2+100" "$1":"$3-100"-"$3 "| grep -v \> ")}}' \
            > {output}
        """


rule ciri2_ref_bsj_blastn:
    input:
        ref = os.path.join(dirs["results"],"ciri2","ref_bsj_seqs.fa"),
        bsj = os.path.join(dirs["results"],"ciri2","bsj_assembly.fasta")
    output:
        os.path.join(dirs["results"], "ciri2", "ref_bsj_blastn.outfmt6")
    conda:
        os.path.join(dirs["envs"], "blast.yaml")
    shell:
        """
        blastn -query {input.bsj} -subject {input.ref} -outfmt 6 \
            | sort -k1,1 -k2,2 -k7,7n \
            > {output}
        """


rule ciri2_bsj_coords:
    input:
        os.path.join(dirs["results"],"ciri2","ref_bsj_blastn.outfmt6")
    output:
        os.path.join(dirs["results"],"ciri2","bsj_coords.tsv")
    log:
        os.path.join(dirs["logs"], "ciri2_bsj_coords.err")
    script:
        os.path.join(dirs["scripts"], "ciri2_blast_to_coords.py")


rule ciri2_bsj_junction_seqs:
    input:
        bsj = os.path.join(dirs["results"],"ciri2","bsj_assembly.fasta"),
        coords = os.path.join(dirs["results"],"ciri2","bsj_coords.tsv")
    output:
        os.path.join(dirs["results"],"ciri2","bsj_junction_seqs.fasta")
    conda:
        os.path.join(dirs["envs"],"star.yaml")
    shell:
        """
        cat {input.coords} \
            | awk '{{print ">"$3"_"$1"|"$2; system("samtools faidx {input.bsj} " $1":"$2-50"-"$2+50" | grep -v \> ")}}' \
            > {output}
        """