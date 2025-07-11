
import logging
import os
import sys
import re
import pandas as pd
import numpy as np
import gzip
import pickle


def detect_bsj_blastn_outfmt6(blast_outfmt6, min_align_len=50, max_bsj_gap=5, min_pident=95, valid_coord_dist=10):
    blast_cols = [
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
    ]
    blast_df = pd.read_csv(blast_outfmt6, sep='\t', header=None, names=blast_cols)
    blast_df = blast_df[
        (blast_df['length'] >= min_align_len) &
        (blast_df['pident'] >= min_pident)
        ]

    blast_df = blast_df.sort_values(by=['qseqid', 'qstart']).reset_index(drop=True)

    bsj_hits = dict()
    val_hits = dict()

    for q_id, group in blast_df.groupby('qseqid'):
        # Further group by subject (chromosome/contig)
        if not q_id in bsj_hits:
            bsj_hits[q_id] = []
        if not q_id in val_hits:
            val_hits[q_id] = []
        for s_id, sub_group in group.groupby('sseqid'):
            # We need at least two alignments to the same reference chromosome/contig to form a BSJ
            if len(sub_group) < 2:
                continue
            # Iterate through all unique pairs of alignments within the same query and subject
            for i in range(len(sub_group)):
                for j in range(i + 1, len(sub_group)):
                    aln1 = sub_group.iloc[i]
                    aln2 = sub_group.iloc[j]

                    # Normalize query coordinates (qstart always < qend)
                    q1_coord_start, q1_coord_end = sorted([aln1['qstart'], aln1['qend']])
                    q2_coord_start, q2_coord_end = sorted([aln2['qstart'], aln2['qend']])

                    # skip large gap or overlap on assembled seq
                    if abs(q2_coord_start - q1_coord_end) > max_bsj_gap:
                        continue

                    # Normalize subject coordinates (sstart always < send for range)
                    s1_coord_start, s1_coord_end = sorted([aln1['sstart'], aln1['send']])
                    s2_coord_start, s2_coord_end = sorted([aln2['sstart'], aln2['send']])

                    # skip large gap or overlap on ref seq
                    if abs(s2_coord_start - s1_coord_end) > max_bsj_gap:
                        continue

                    # raw query coordinates
                    q1_start, q1_end = [aln1['qstart'], aln1['qend']]
                    q2_start, q2_end = [aln2['qstart'], aln2['qend']]

                    # raw subject coordinates
                    s1_start, s1_end = [aln1['sstart'], aln1['send']]
                    s2_start, s2_end = [aln2['sstart'], aln2['send']]

                    # Reverse alignment
                    if s1_start > s1_end and s2_start > s2_end:
                        # sanity check query direction
                        if q1_start < q1_end and q2_start < q2_end:
                            # bsj condition
                            if q2_start > q1_start and s2_start > s1_start:
                                bsj_hits[q_id].append([q1_end, s_id])
                            else:
                                val_hits[q_id].append([q1_end, s_id])
                    # Forward alignment
                    elif s1_start < s1_end and s2_start < s1_end:
                        # sanity check query direction
                        if q1_start < q1_end and q2_start < q2_end:
                            # bsj condition
                            if q2_start > q1_start and s2_start > s1_start:
                                bsj_hits[q_id].append([q1_end, s_id])
                            else:
                                val_hits[q_id].append([q1_end, s_id])
        for bsj_hit_coord in bsj_hits[q_id]:
            is_bsj = True
            for val_hit_coord in val_hits[q_id]:
                if abs(bsj_hit_coord[0] - val_hit_coord[0]) <= valid_coord_dist:
                    is_bsj = False
                    break
            if is_bsj:
                yield f"{q_id}\t{bsj_hit_coord[0]}\t{bsj_hit_coord[1]}\n"



def main(snakemake):

    logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)

    unique_bsj_hits = set()
    for bsj_hit in detect_bsj_blastn_outfmt6(snakemake.input[0]):
        unique_bsj_hits.add(bsj_hit)
    with open(snakemake.output[0], "w") as outfh:
        for bsj_hit in unique_bsj_hits:
            outfh.write(bsj_hit)



if __name__ == "__main__":
    main(snakemake)