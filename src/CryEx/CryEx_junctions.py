#!/usr/bin/env python3
import logging
import multiprocessing as mp
import optparse
import os
import subprocess
from collections import defaultdict

import pandas as pd
import pysam


JunctionKey = tuple[str, int, int, str]


def read_filtering(read: pysam.AlignedSegment) -> bool:
    return not (
        read.is_unmapped
        or read.is_secondary
        or read.is_supplementary
        or read.is_duplicate
        or read.is_qcfail
    )


def determine_strand(read: pysam.AlignedSegment, strandedness: str) -> str:
    strand_pairs = {
        "fr-firststrand": ("-", "+"),
        "fr-secondstrand": ("+", "-"),
        "unstranded": (".", "."),
    }
    try:
        strand1, strand2 = strand_pairs[strandedness]
    except KeyError as error:
        raise ValueError("Strandedness is not valid") from error

    if read.is_read1 and not read.is_reverse:
        return strand1
    if read.is_read1 and read.is_reverse:
        return strand2
    if read.is_read2 and read.is_reverse:
        return strand1
    if read.is_read2 and not read.is_reverse:
        return strand2
    return strand1 if not read.is_reverse else strand2


def extract_jxns_from_bam(
    bam_file: str,
    star_file: str,
    sample_name: str,
    strandedness: str,
    overhang_minimum: int = 8,
) -> tuple[str, dict[JunctionKey, int]]:
    junctions: defaultdict[JunctionKey, int] = defaultdict(int)

    if not os.path.exists(star_file):
        logging.info(
            "STAR junction file %s does not exist; counting junctions in %s",
            star_file,
            bam_file,
        )
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for read in bam:
                if not read_filtering(read):
                    continue

                read_strand = determine_strand(read, strandedness)
                blocks = read.get_blocks()
                chrom = read.reference_name
                if chrom is None:
                    continue

                for index in range(len(blocks) - 1):
                    junction_start = blocks[index][1]
                    junction_end = blocks[index + 1][0]
                    if junction_end - junction_start > 20:
                        key = (chrom, junction_start, junction_end, read_strand)
                        junctions[key] += 1
    else:
        dataframe = pd.read_csv(star_file, header=None, sep="\t")
        dataframe.columns = [
            "chrom",
            "jxn_s",
            "jxn_e",
            "strand",
            "motif",
            "annotated",
            "n_uniq",
            "n_multi",
            "max_ov",
        ]
        dataframe = dataframe[dataframe["max_ov"] >= overhang_minimum]
        strand_lookup = {0: ".", 1: "+", 2: "-"}

        for row in dataframe.itertuples():
            if row.strand not in strand_lookup:
                raise ValueError(f"Strand not valid: {row.strand}")
            junction = (
                row.chrom,
                row.jxn_s - 1,
                row.jxn_e,
                strand_lookup[row.strand],
            )
            junctions[junction] = row.n_uniq

    return sample_name, dict(junctions)


def main() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(filename)s %(asctime)s %(levelname)s:\t%(message)s",
    )
    usage = "CryEx_junctions -f fofn -o output_file"
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

    logging.info("BAM file strandedness = %s", options.strandedness)
    metadata = pd.read_csv(options.fofn, sep="\t")
    task_args = [
        (
            str(row.BAM),
            str(row.STAR_SJ_OUT),
            str(row.SAMPLE),
            options.strandedness,
        )
        for row in metadata.itertuples()
    ]

    with mp.Pool(processes=10) as pool:
        results = pool.starmap(extract_jxns_from_bam, task_args)

    merged: defaultdict[JunctionKey, defaultdict[str, int]] = defaultdict(
        lambda: defaultdict(int)
    )
    for sample_name, sample_junctions in results:
        for junction, count in sample_junctions.items():
            merged[junction][sample_name] = count

    sample_names = metadata["SAMPLE"].astype(str).tolist()
    with open(options.output, "w") as output:
        for (chrom, start, end, strand), sample_counts in merged.items():
            counts = [sample_counts[sample] for sample in sample_names]
            output.write(
                f"{chrom}\t{start}\t{end}\t{','.join(sample_names)}\t"
                f"{','.join(map(str, counts))}\t{strand}\n"
            )

    # Sort and compress without a shell so paths are passed literally.
    with open(options.output, "rb") as junction_input:
        sort_process = subprocess.Popen(
            ["sort", "-k1,1V", "-k2,2n", "-k3,3n"],
            stdin=junction_input,
            stdout=subprocess.PIPE,
        )

    try:
        with open(f"{options.output}.gz", "wb") as compressed_output:
            assert sort_process.stdout is not None
            bgzip_process = subprocess.Popen(
                ["bgzip", "-c"],
                stdin=sort_process.stdout,
                stdout=compressed_output,
            )
            sort_process.stdout.close()
            bgzip_returncode = bgzip_process.wait()
        sort_returncode = sort_process.wait()
    except Exception:
        if sort_process.stdout is not None:
            sort_process.stdout.close()
        sort_process.terminate()
        sort_process.wait()
        raise

    if sort_returncode != 0:
        raise subprocess.CalledProcessError(sort_returncode, "sort")
    if bgzip_returncode != 0:
        raise subprocess.CalledProcessError(bgzip_returncode, "bgzip")

    subprocess.run(
        ["tabix", "-p", "bed", "-f", f"{options.output}.gz"], check=True
    )
    logging.info("Junction processing finished")


if __name__ == "__main__":
    main()
