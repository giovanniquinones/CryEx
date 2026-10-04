#!/usr/bin/env python3
import logging
import optparse
import subprocess

import pandas as pd


def main() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(filename)s %(asctime)s %(levelname)s:\t%(message)s",
    )
    usage = "CryEx_stringtie -f fofn -o output_file"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-f", dest="fofn", help="Metadata file [required]")
    parser.add_option(
        "-s",
        dest="strandedness",
        help="Library strandedness [optional]",
        default="fr-firststrand",
    )
    parser.add_option("-o", dest="output", help="Output file [required]")
    options, _ = parser.parse_args()

    if not options.fofn or not options.output:
        parser.error(f"Usage: {usage}")

    strand_params = {
        "fr-firststrand": ["--rf"],
        "fr-secondstrand": ["--fr"],
        "unstranded": [],
    }
    if options.strandedness not in strand_params:
        parser.error(
            "Strandedness must be fr-firststrand, fr-secondstrand, or unstranded"
        )

    metadata = pd.read_csv(options.fofn, sep="\t")
    sample_gtf_files: list[str] = []
    for row in metadata.itertuples():
        sample_gtf = f"{options.output}.{row.SAMPLE}.gtf"
        sample_gtf_files.append(sample_gtf)

        try:
            with open(sample_gtf, "w") as output_handle:
                result = subprocess.run(
                    [
                        "stringtie",
                        str(row.BAM),
                        "-j",
                        "2",
                        "-f",
                        "0.05",
                        "-g",
                        "10",
                        *strand_params[options.strandedness],
                    ],
                    stdout=output_handle,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True,
                )
        except subprocess.CalledProcessError as error:
            logging.error("StringTie failed for sample %s", row.SAMPLE)
            logging.error("Return code: %s", error.returncode)
            logging.error("Error message: %s", error.stderr)
            raise

        if result.stderr:
            logging.info(result.stderr.rstrip())

    with open(options.output, "w") as output_handle:
        result = subprocess.run(
            ["sort", "-k1,1V", "-k4,4n", "-k5,5n", "--", *sample_gtf_files],
            stdout=output_handle,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )
    if result.stderr:
        logging.info(result.stderr.rstrip())
    logging.info("Combined GTF written to %s", options.output)


if __name__ == "__main__":
    main()
