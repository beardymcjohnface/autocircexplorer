#!/usr/bin/env python3


"""Convert GTF for flatfile format


"""

import logging
import pandas as pd


def parse_salmon_quants(resdir, samples, filename):
    for sample in samples:
        df = pd.read_table(os.path.join(resdir, sample, filename))
        df = df[["Name", "TPM", "NumReads"]]
        df["Sample"] = sample
        yield df


def main(snakemake):
    logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)

    logging.debug("Slurping Salmon outputs")
    all_sample_df = pd.concat(
        parse_salmon_quants(snakemake.params.dir, snakemake.params.samples, snakemake.params.file)
    )

    logging.debug("Pivot for TPMs")
    tpm_wide = all_sample_df.pivot("Sample", index= "Name", values = "TPM")

    logging.debug("Pivot for NumReads")
    numreads_wide = all_sample_df.pivot("Sample", index="Name", values="NumReads")

    logging.debug("Saving wide tables")
    tpm_wide.to_csv(
        snakemake.output.tpm,
        index=False,
        compression={"method": "gzip", "compresslevel": 1, "mtime": 1}
    )
    numreads_wide.to_csv(
        snakemake.output.num,
        index=False,
        compression={"method": "gzip", "compresslevel": 1, "mtime": 1}
    )

    logging.debug("Done!")


if __name__ == "__main__":
    main(snakemake)
