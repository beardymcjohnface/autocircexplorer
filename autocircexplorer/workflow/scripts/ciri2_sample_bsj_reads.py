
import logging
import os
import sys
import re
import pandas as pd
import gzip
import pickle


def parse_ciri2(ciri2_output, min_reads):
    """
    return (dict): read_id = bsj_id
                   use the dict keys to pull reads for all BSJs in one pass
                   and filter into their individual junctions
    """
    out_read_ids = dict()
    with open(ciri2_output, "r") as fh:
        _ = fh.readline()
        for line in fh:
            l = line.split("\t")
            if int(l[4]) < min_reads:
                read_ids = l[11].strip(",\n").split(",")
                for read_id in read_ids:
                    out_read_ids["@" + read_id] = l[0]
    return out_read_ids


def parse_fastq(filepath):
    if filepath.endswith('.gz'):
        _open = gzip.open
        mode = 'rt'
    else:
        _open = open
        mode = 'r'

    with _open(filepath, mode) as f:
        while True:
            header_id = f.readline().strip()
            if not header_id:
                break
            sequence = f.readline().strip()
            quality_header = f.readline().strip()
            quality_scores = f.readline().strip()

            yield header_id, sequence, quality_header, quality_scores


def pull_reads(bsj_ids, fq_file):
    pull_ids = set(bsj_ids.keys())
    out_bsj_reads = dict()
    for fq_id, fq_seq, fq_qh, fq_qs in parse_fastq(fq_file):
        fq_match = fq_id.split(" ")[0]
        if fq_match in pull_ids:
            if not bsj_ids[fq_match] in out_bsj_reads.keys():
                out_bsj_reads[bsj_ids[fq_match]] = list()
            out_bsj_reads[bsj_ids[fq_match]].append(f"{fq_id}\n{fq_seq}\n{fq_qh}\n{fq_qs}\n")
    return out_bsj_reads


def main(snakemake):

    logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)
    bsj_read_ids = parse_ciri2(snakemake.input.ciri2, snakemake.params.min_reads)
    bsj_reads = dict()
    bsj_reads["R1"] = pull_reads(bsj_read_ids, snakemake.input.r1)
    if snakemake.params.r2:
        bsj_reads["R2"] = pull_reads(bsj_read_ids, snakemake.params.r2)

    with open(snakemake.output[0], "wb") as outfh:
        pickle.dump(bsj_reads, outfh, protocol=pickle.HIGHEST_PROTOCOL)



if __name__ == "__main__":
    main(snakemake)