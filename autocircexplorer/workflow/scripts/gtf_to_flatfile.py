#!/usr/bin/env python3


"""Convert GTF for flatfile format


"""

import logging
import re


def initialise_transcript(geneName, chrom, strand):
    return {
        "geneName": geneName,
        "isoformName": "",
        "chrom": chrom,
        "strand": strand,
        "exonStarts": list(),
        "exonEnds": list(),
        "exonCount": 0
    }


def write_transcript(transcript_dict, fh):
    if transcript_dict["isoformName"] != "":
        logging.debug("Writing transcript: " + transcript_dict["isoformName"])

        fh.write("\t".join([
            transcript_dict["geneName"],
            transcript_dict["isoformName"],
            transcript_dict["chrom"],
            transcript_dict["strand"],
            transcript_dict["txStart"],
            transcript_dict["txEnd"],
            transcript_dict["cdsStart"],
            transcript_dict["cdsEnd"],
            str(transcript_dict["exonCount"]),
            ",".join(transcript_dict["exonStarts"]),
            ",".join(transcript_dict["exonEnds"]) + "\n",
        ]))


def new_transcript(transcript_dict):
    transcript_dict["exonCount"] = 0
    transcript_dict["exonStarts"] = list()
    transcript_dict["exonEnds"] = list()
    return transcript_dict


def parse_gtf(snakemake):
    logging.debug("Reading GTF")
    with open(snakemake.input[0], "r") as in_gtf:

        logging.debug("opening output for writing")
        with open(snakemake.output[0], "w") as out_flat:
            transcript_dict = initialise_transcript("","", "")

            for line in in_gtf:
                fields = line.split("\t")

                if len(fields) < 9:
                    continue

                if fields[2] == "gene":
                    write_transcript(transcript_dict, out_flat)
                    geneName = re.search(r'; gene_name "([^"]*)";', fields[8]).group(1)
                    transcript_dict = initialise_transcript(geneName, fields[0], fields[6])
                    logging.debug("New gene: " + geneName)

                elif fields[2] == "transcript":
                    write_transcript(transcript_dict, out_flat)
                    transcript_dict = new_transcript(transcript_dict)
                    isoformName = re.search(r'; transcript_id "([^"]*)";', fields[8]).group(1)
                    transcript_dict["isoformName"] = isoformName
                    transcript_dict["txStart"] = fields[3]
                    transcript_dict["txEnd"] = fields[4]
                    transcript_dict["cdsStart"] = fields[3]
                    transcript_dict["cdsEnd"] = fields[4]

                elif fields[2] == "exon":
                    transcript_dict["exonStarts"].append(fields[3])
                    transcript_dict["exonEnds"].append(fields[4])

            write_transcript(transcript_dict, out_flat)

def main(snakemake):
    logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)

    logging.debug("Parsing GTF and printing flat file on the fly")
    parse_gtf(snakemake)

    logging.debug("Done!")



if __name__ == "__main__":
    main(snakemake)
