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


def collect_rmats_multi_out(snakemake):
    libcounts = pd.read_table(snakemake.input.lib, header=None, index_col=0)
    libcounts = libcounts.to_dict()[1]

    for i in ["A3SS", "A5SS", "MXE", "RI", "SE"]:
        df = pd.read_table(os.path.join(snakemake.params.outpath, i + ".MATS." + snakemake.params.type + ".txt"))

        df[[ "S1:IJC:" + x for x in snakemake.params.s1 ]] = df["IJC_SAMPLE_1"].str.split(",", expand=True)
        df[[ "S2:IJC:" + x for x in snakemake.params.s2 ]] = df["IJC_SAMPLE_2"].str.split(",", expand=True)
        df[[ "S1:SJC:" + x for x in snakemake.params.s1 ]] = df["SJC_SAMPLE_1"].str.split(",", expand=True)
        df[[ "S2:SJC:" + x for x in snakemake.params.s2 ]] = df["SJC_SAMPLE_2"].str.split(",", expand=True)

        df.drop(columns=["IJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_1", "SJC_SAMPLE_2"], inplace=True)

        df.to_csv(os.path.join(snakemake.params.respath, i + ".raw.tsv"), sep='\t', index=False, float_format='%.6f')

        cols_for_long_output = list()

        for sample in snakemake.params.s1:
            cols_for_long_output.append("S1:IJC:" + sample)
            cols_for_long_output.append("S1:SJC:" + sample)
            df["S1:IJC:" + sample] = (pd.to_numeric(df["S1:IJC:" + sample]) / libcounts[sample]) * 1000000
            df["S1:SJC:" + sample] = (pd.to_numeric(df["S1:SJC:" + sample]) / libcounts[sample]) * 1000000

        for sample in snakemake.params.s2:
            cols_for_long_output.append("S2:IJC:" + sample)
            cols_for_long_output.append("S2:SJC:" + sample)
            df["S2:IJC:" + sample] = (pd.to_numeric(df["S2:IJC:" + sample]) / libcounts[sample]) * 1000000
            df["S2:SJC:" + sample] = (pd.to_numeric(df["S2:SJC:" + sample]) / libcounts[sample]) * 1000000

        df.to_csv(os.path.join(snakemake.params.respath, i + ".CPM.tsv"), sep='\t', index=False, float_format='%.6f')

        for col_name in cols_for_long_output:
            df_out = df[["ID", col_name]]
            df_out.rename(columns={col_name: "CPM"}, inplace=True)
            sample_group, sample_type, sample_name = col_name.split(":")
            df_out["ID"] = df_out["ID"].astype(str) + ":" + sample_type
            df_out["Sample"] = sample_name
            df_out["Group"] = sample_group
            df_out["Type"] = sample_type
            df_out["Event"] = i
            yield df_out[["ID", "Sample", "Group", "Type", "Event", "CPM"]]


def event_wide_counts(df, event_type):
    df = df.loc[df["eventType"] == event_type]
    df = df[["exonID", "sampleID", "skipCounts", "inclusionCounts"]].groupby(["exonID","sampleID"]).sum().reset_index()
    df_skip = df.pivot_table(index="exonID", columns="sampleID", values="skipCounts", fill_value=0).reset_index()
    df_incl = df.pivot_table(index="exonID", columns="sampleID", values="inclusionCounts", fill_value=0).reset_index()
    return df_skip, df_incl


def main(snakemake):

    logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)

    if snakemake.params.mode == "single":
        logging.debug("Reading in RMATS output files")
        all_outs = collect_rmats_out(snakemake.input.rmats_out)

        logging.debug("Combining RMATS outputs")
        combined_all_outs = pd.concat(all_outs)

        logging.debug("Writing long output")
        combined_all_outs.to_csv(snakemake.output.long, index=False, compression={"method": "gzip", "compresslevel": 1, "mtime": 1})

        logging.debug("Writing wide table skip counts")
        for event_type in ["A3SS", "A5SS", "MXE", "RI", "SE"]:
            event_skip, event_incl = event_wide_counts(combined_all_outs, event_type)
            event_skip.to_csv(
                os.path.join(snakemake.output.outpath, event_type + ".skip.csv.gz"),
                index=False, compression={"method": "gzip", "compresslevel": 1, "mtime": 1}
            )
            event_incl.to_csv(
                os.path.join(snakemake.output.outpath, event_type + ".inclusion.csv.gz"),
                index=False, compression={"method": "gzip", "compresslevel": 1, "mtime": 1}
            )

        logging.debug("Writing exon meta table")
        exon_meta = combined_all_outs[
            ["exonID", "GeneID", "geneSymbol", "chr", "strand", "exonStart", "exonEnd", "eventType"]].drop_duplicates()
        exon_meta.to_csv(
            os.path.join(snakemake.output.outpath, "exon_id_metadata.csv.gz"),
            index=False, compression={"method": "gzip", "compresslevel": 1, "mtime": 1}
        )

    else:
        logging.debug("Reading in RMATS output files")
        all_outs = collect_rmats_multi_out(snakemake)

        logging.debug("Combining RMATS outputs")
        combined_all_outs = pd.concat(all_outs)

        logging.debug("Writing long output")
        combined_all_outs.to_csv(snakemake.output.long, index=False, compression={"method": "gzip", "compresslevel": 1, "mtime": 1})



if __name__ == "__main__":
    main(snakemake
        # snakemake.input.rmats_out,
        # snakemake.output.long,
        # snakemake.params.outpath,
        # snakemake.params.mode,
        # snakemake.log[0],
    )
