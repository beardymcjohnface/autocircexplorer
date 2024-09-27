import glob
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
    samples["group2"] = list(samples["reads"].keys())
    samples["reads"].update(fastq_finder.parse_samples_to_dictionary(config["args"]["input1"]))
    samples["group1"] = list(set(samples["reads"].keys()) - set(samples["group2"]))
else:
    samples["reads"] = fastq_finder.parse_samples_to_dictionary(config["args"]["input1"])
    samples["group1"] = []
    samples["group2"] = []

samples["names"] = list(samples["reads"].keys())


# Targets
targets = [
    expand(
        os.path.join(dirs["results"], "star", "{sample}.bam"),
        sample=samples["names"]
    ),
    expand(
        os.path.join(dirs["results"], "ce2", "{sample}.circexplorer2.parse"),
        sample=samples["names"]
    )
]

if config["args"]["input2"]:
    targets.append(os.path.join(dirs["results"], "rmats_multi_summary_long.tsv"))
else:
    targets.append(
        expand(
            os.path.join(dirs["results"], "{event}.{count}.csv.gz"),
            event=["A3SS","A5SS","MXE","RI","SE"],
            count=["skip","inclusion"]
        ),
    )


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

    files = glob.glob(os.path.join(".snakemake", "log", "*.snakemake.log"))
    if files:
        current_log = max(files, key=os.path.getmtime)
        shell("cat " + current_log + " >> " + config["args"]["log"])


onsuccess:
    copy_log_file()

onerror:
    copy_log_file()
