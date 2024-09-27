# autoCircExplorer

Run CIRCexplorer2 and rmats-turbo but better.

## Install

```shell
# OPTIONAL: set up an env
conda create -n autoce2 python=3.12
conda activate autoce2

# clone the repo
git clone https://github.com/beardymcjohnface/autocircexplorer.git

# install with pip
cd autocircexplorer/
pip install -e .
```

## Run

__Run rmats individually__

Make sure your sequencing reads are all in one directory, or pass a sample.tsv file.

```shell
autoce2 run \
    --input1 fastqDir/ \
    --ref starRefDir/ \
    --gtf refGenes.gtf
```

__Compare two groups__

Pass separate directories or sample files for each sample group.

```shell
autoce2 run \
    --input1 group1Fastq/ \
    --input2 group2Fastq/ \
    --ref starRefDir/ \
    --gtf refGenes.gtf
```
