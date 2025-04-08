import os

from metasnek import fastq_finder, fasta_finder


# Directories
dirs = {
    "logs": os.path.join(config["args"]["output"], "logs"),
    "bench": os.path.join(config["args"]["output"], "bench"),
    "results": os.path.join(config["args"]["output"], "results"),
    "envs": os.path.join(workflow.basedir, "envs"),
    "scripts": os.path.join(workflow.basedir, "scripts")
}


# PARSE SAMPLES
samples = dict()

if config["args"]["input2"]:
    samples["reads"] = fastq_finder.parse_samples_to_dictionary(config["args"]["input2"])
    samples["group2"] = sorted(list(samples["reads"].keys()))
    samples["reads"].update(fastq_finder.parse_samples_to_dictionary(config["args"]["input1"]))
    samples["group1"] = sorted(list(set(samples["reads"].keys()) - set(samples["group2"])))
else:
    samples["reads"] = fastq_finder.parse_samples_to_dictionary(config["args"]["input1"])
    samples["group1"] = []
    samples["group2"] = []

samples["names"] = list(samples["reads"].keys())


# Targets
targets = dict()
targets["star"] = [
    expand(
        os.path.join(dirs["results"], "star", "{sample}.{bam}"),
        sample=samples["names"],
        bam=["bam","bam.bai"]
    )
]

targets["flagstat"] = [
    expand(
        os.path.join(dirs["results"], "multiqc", "{sample}.Log.final.out"),
        sample=samples["names"]
    )
]

targets["salmon"] = [
    expand(
        os.path.join(dirs["results"], "salmon", "{sample}", "quant.sf"),
        sample=samples["names"]
    ),
    os.path.join(dirs["results"],"salmon","TPM.tsv"),
    os.path.join(dirs["results"],"salmon","NumReads.tsv")
]

targets["ce2"] = [
    expand(
        os.path.join(dirs["results"], "ce2", "{sample}.circexplorer2.{step}"),
        sample=samples["names"],
        step = ["parse", "annotated"]
    ),
    os.path.join(dirs["results"], "lib.counts.tsv"),
    os.path.join(dirs["results"], "ce2", "refFlatFile.txt")
]

targets["ciri2"] = expand(
    os.path.join(dirs["results"], "ciri2", "{sample}.ciri2"),
    sample=samples["names"]
)

if config["args"]["input2"]:
    targets["rmats"] = [
        expand(
            os.path.join(dirs["results"], "{file}.{count}.tsv"),
            file=["A3SS", "A5SS", "MXE", "RI", "SE"],
            count=["raw","CPM"]
        ),
        os.path.join(dirs["results"],"rmats_multi_summary_long.csv.gz")
    ]
else:
    targets["rmats"] = expand(
        os.path.join(dirs["results"], "{event}.{count}.csv.gz"),
        event=["A3SS","A5SS","MXE","RI","SE"],
        count=["skip","inclusion"]
    )

targets["fastqc"] = list()
for sample in samples["names"]:
    targets["fastqc"].append(os.path.join(dirs["results"], "multiqc", sample + "_R1_fastqc.zip"))
    if samples["reads"][sample]["R2"]:
        targets["fastqc"].append(os.path.join(dirs["results"],"multiqc",sample + "_R2_fastqc.zip"))


# Misc
target_rules = []



def targetRule(fn):
    """Mark rules as target rules for rule print_targets"""
    assert fn.__name__.startswith("__")
    target_rules.append(fn.__name__[2:])
    return fn


def copy_log_file():
    """Concatenate Snakemake log to output log file"""
    import glob

    def rename_file_if_exists(basename, counter):
        if counter==0:
            exist_file = basename
        else:
            exist_file = basename + '.' + str(counter)
        if os.path.exists(exist_file):
            rename_file_if_exists(basename,counter + 1)
            os.rename(exist_file,basename + '.' + str(counter + 1))

    if os.path.exists(config["args"]["log"]):
        rename_file_if_exists(config["args"]["log"], 0)
    files = glob.glob(os.path.join(".snakemake", "log", "*.snakemake.log"))
    if files:
        current_log = max(files, key=os.path.getmtime)
        shell("cat " + current_log + " > " + config["args"]["log"])


onsuccess:
    copy_log_file()

onerror:
    copy_log_file()
