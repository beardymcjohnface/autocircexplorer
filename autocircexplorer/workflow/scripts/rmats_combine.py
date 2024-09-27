#!/usr/bin/env python3


"""Combine the output of rmats-turbo


"""


import logging
import os
import sys
import re
import pandas as pd


def rmats_output_specific_cols(in_file_path):
    keep_cols = {
        "A3SS": ["shortES", "shortEE"],
        "A5SS": ["shortES", "shortEE"],
        "SE": ["exonStart_0base", "exonEnd"],
        "MXE": ["1stExonStart_0base", "1stExonEnd"],
        "RI": ["riExonStart_0base", "riExonEnd"]
    }
    rmats_type = os.path.basename(in_file_path).replace(".MATS.JC.txt", "").replace(".MATS.JCEC.txt", "")
    return rmats_type, keep_cols[rmats_type]


def convert_df(df):
    df = df.rename(columns = {
        "shortES": "exonStart",
        "shortEE": "exonEnd",
        "exonStart_0base": "exonStart",
        "1stExonStart_0base": "exonStart",
        "1stExonEnd": "exonEnd",
        "riExonStart_0base": "exonStart",
        "riExonEnd": "exonEnd",
         "IJC_SAMPLE_1": "inclusionCounts",
         "SJC_SAMPLE_1": "skipCounts"}
    )
    df["exonStart"] = df["exonStart"].astype(str)
    df["exonEnd"] = df["exonEnd"].astype(str)
    df["exonID"] = df["chr"] + ":" + df["exonStart"] + "-" + df["exonEnd"] + ":" + df["eventType"]
    return df


def collect_rmats_out(input_files):

    common_cols_to_keep = ["GeneID", "geneSymbol", "chr", "strand", "IJC_SAMPLE_1", "SJC_SAMPLE_1"]

    for rmats_file in input_files:
        rmats_type, specific_cols_to_keep = rmats_output_specific_cols(rmats_file)
        cols_to_keep = common_cols_to_keep + specific_cols_to_keep
        rmats_df = pd.read_table(rmats_file, sep="\t", usecols=cols_to_keep)
        rmats_df["eventType"] = rmats_type
        rmats_df = convert_df(rmats_df)

        sample_id  = re.search(r'rmats/(.*).rmats_single', rmats_file).group(1)
        rmats_df["sampleID"] = sample_id
        yield rmats_df


def event_wide_counts(df, event_type):
    df = df.loc[df["eventType"] == event_type]
    df = df[["exonID", "sampleID", "skipCounts", "inclusionCounts"]].groupby(["exonID","sampleID"]).sum().reset_index()
    df_skip = df.pivot_table(index="sampleID", columns="exonID", values="skipCounts", fill_value=0).reset_index()
    df_incl = df.pivot_table(index="sampleID", columns="exonID", values="inclusionCounts", fill_value=0).reset_index()
    return df_skip, df_incl


def main(input_files, out_long, out_path, log_file, **kwargs):

    logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG)

    logging.debug("Reading in RMATS output files")
    all_outs = collect_rmats_out(input_files)

    logging.debug("Combining RMATS outputs")
    combined_all_outs = pd.concat(all_outs)

    logging.debug("Writing long output")
    combined_all_outs.to_csv(out_long, index=False, compression={"method": "gzip", "compresslevel": 1, "mtime": 1})

    logging.debug("Writing wide table skip counts")
    for event_type in ["A3SS","A5SS","MXE","RI","SE"]:
        event_skip, event_incl = event_wide_counts(combined_all_outs, event_type)
        event_skip.to_csv(
            os.path.join(out_path, event_type + ".skip.csv.gz"),
            index=False, compression={"method": "gzip", "compresslevel": 1, "mtime": 1}
        )
        event_incl.to_csv(
            os.path.join(out_path, event_type + ".inclusion.csv.gz"),
            index=False, compression={"method": "gzip", "compresslevel": 1, "mtime": 1}
        )

    logging.debug("Writing exon meta table")
    exon_meta = combined_all_outs[["exonID", "GeneID", "geneSymbol", "chr", "strand", "exonStart", "exonEnd", "eventType"]].drop_duplicates()
    exon_meta.to_csv(
        os.path.join(out_path, "exon_id_metadata.csv.gz"),
        index=False, compression={"method": "gzip", "compresslevel": 1, "mtime": 1}
    )


if __name__ == "__main__":
    main(
        snakemake.input,
        snakemake.output.long,
        snakemake.params.outpath,
        snakemake.log[0],
    )
