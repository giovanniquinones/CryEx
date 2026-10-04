#!/usr/bin/env python3
import logging
import optparse
from collections import Counter

import numpy as np
import pandas as pd
import pysam


def read_filtering(pysam_read: pysam.AlignedSegment) -> bool:
    if pysam_read.is_unmapped:
        return False
    if pysam_read.is_secondary:
        return False
    if pysam_read.is_supplementary:
        return False
    if pysam_read.is_duplicate:
        return False
    if pysam_read.is_qcfail:
        return False
    return True


def determine_strand(read: pysam.AlignedSegment, strandedness: str) -> str:
    if strandedness == "fr-firststrand":
        strand1, strand2 = "-", "+"
    elif strandedness == "fr-secondstrand":
        strand1, strand2 = "+", "-"
    elif strandedness == "unstranded":
        strand1, strand2 = ".", "."
    else:
        raise ValueError("Strandedness is not valid")

    if read.is_read1 and not read.is_reverse:
        return strand1
    if read.is_read1 and read.is_reverse:
        return strand2
    if read.is_read2 and read.is_reverse:
        return strand1
    if read.is_read2 and not read.is_reverse:
        return strand2
    return strand1 if not read.is_reverse else strand2


def splicing_type(
    jun5ss_counter: Counter[int], jun3ss_counter: Counter[int], strand: str
) -> str:
    has_5ss = bool(jun5ss_counter)
    has_3ss = bool(jun3ss_counter)

    if not has_5ss and not has_3ss:
        return "not_spliced"
    if has_5ss and has_3ss:
        return "middle_exon"
    if strand == ".":
        return "end_exon"
    if (has_5ss and strand == "-") or (has_3ss and strand == "+"):
        return "last_exon"
    if (has_5ss and strand == "+") or (has_3ss and strand == "-"):
        return "first_exon"
    return "unknown"


def get_coverage(
    chrom: str,
    exon_start: int,
    exon_end: int,
    strand: str,
    bh: pysam.AlignmentFile,
    strandedness: str,
    min_overlap: int = 5,
    padding: int = 10,
) -> tuple[int, str, float, int]:
    read_lengths: Counter[int] = Counter()
    jun5ss_counter: Counter[int] = Counter()
    jun3ss_counter: Counter[int] = Counter()
    inclusion_n = 0
    exon_start0 = exon_start - 1

    for read in bh.fetch(chrom, exon_start0, exon_end + 1):
        if not read_filtering(read):
            continue
        if read.get_overlap(exon_start0, exon_end) < min_overlap:
            continue

        read_strand = determine_strand(read, strandedness)
        if strandedness != "unstranded" and read_strand != strand:
            continue

        inclusion_n += 1
        blocks = read.get_blocks()
        for index in range(len(blocks) - 1):
            jun5ss = blocks[index][1]
            jun3ss = blocks[index + 1][0]

            if exon_end - padding <= jun5ss <= exon_end + padding:
                jun5ss_counter[jun5ss] += 1
            if exon_start0 - padding <= jun3ss <= exon_start0 + padding:
                jun3ss_counter[jun3ss] += 1

        if read.query_length is not None:
            read_lengths[read.query_length] += 1

    read_length = read_lengths.most_common(1)[0][0] if read_lengths else np.nan
    label = splicing_type(jun5ss_counter, jun3ss_counter, strand)
    max_junc_count = max([*jun5ss_counter.values(), *jun3ss_counter.values(), 0])
    return inclusion_n, label, read_length, max_junc_count


def calculate_psi(
    inclusion_n: int, exclusion_n: int, interval_length: int, read_length: float
) -> float:
    intron_norm = inclusion_n / (interval_length + read_length)
    intron_norm -= intron_norm % 0.01

    ex_norm = exclusion_n / read_length
    ex_norm -= ex_norm % 0.01

    try:
        psi = round(intron_norm / (intron_norm + ex_norm), 2)
    except (ZeroDivisionError, TypeError):
        psi = np.nan

    return psi


class JunctionRecord:
    """Parsed junction row from the tabix-indexed junction file."""

    def __init__(self, line: str) -> None:
        (
            contig,
            junction_start,
            junction_end,
            samples,
            counts,
            strand,
        ) = line.split("\t")
        self.start: int = int(junction_start)
        self.end: int = int(junction_end)
        self.contig: str = contig
        self.strand: str = {"1": "+", "2": "-", "0": "."}.get(strand, strand)
        self.length: int = self.end - self.start + 1
        self.samples: list[str] = samples.split(",")
        self.counts: list[int] = list(map(int, counts.split(",")))
        self.sample_counts: dict[str, int] = dict(zip(self.samples, self.counts))


def main() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(filename)s %(asctime)s %(levelname)s:\t%(message)s",
    )
    usage_instructions = (
        "CryEx_psi_calculator -f fofn -e exon_coordinates -o output_file"
    )
    parser = optparse.OptionParser(usage=usage_instructions)
    parser.add_option("-f", dest="fofn", help="Metadata file [required]")
    parser.add_option("-j", dest="jxn_file", help="Junction coordinates file")
    parser.add_option("-e", dest="exons_file", help="Exon coordinates file")
    parser.add_option(
        "-s",
        dest="strandedness",
        help="Library strandedness [optional]",
        default="fr-firststrand",
    )
    parser.add_option("-o", dest="output_file", help="Output file")
    options, _ = parser.parse_args()

    if not all(
        (options.fofn, options.jxn_file, options.exons_file, options.output_file)
    ):
        parser.error(f"Usage: {usage_instructions}")

    metadata = pd.read_csv(options.fofn, sep="\t")
    bam_list = [
        (row.SAMPLE, pysam.AlignmentFile(row.BAM, "rb"))
        for row in metadata.itertuples()
    ]
    exons = pd.read_csv(options.exons_file, sep="\t", header=None, comment="#")
    exons.columns = [
        "chrom",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attributes",
    ]
    exons = exons.loc[
        exons["type"] == "exon", ["chrom", "start", "end", "strand"]
    ].drop_duplicates()

    with open(options.output_file, "w") as output:
        output.write(
            "exon_type\tchrom\texon_3ss\texon_5ss\tstrand\tinclusion_n\t"
            "jxn_inc_n\texc_5ss\texc_3ss\texclusion_n\tSAMPLE\tPSI\n"
        )

        logging.info("Calculating PSI for exons")
        with pysam.TabixFile(options.jxn_file) as junctions:
            for exon_index, exon in enumerate(exons.itertuples()):
                coverage = {
                    sample_name: get_coverage(
                        exon.chrom,
                        exon.start,
                        exon.end,
                        exon.strand,
                        bam,
                        options.strandedness,
                    )
                    for sample_name, bam in bam_list
                }
                exon_length = exon.end - exon.start + 1

                for junction_text in junctions.fetch(
                    exon.chrom, exon.start, exon.end
                ):
                    junction = JunctionRecord(junction_text)
                    if junction.strand != exon.strand:
                        continue
                    if not (
                        junction.start < exon.start < exon.end < junction.end
                    ):
                        continue

                    for sample_name, exclusion_n in junction.sample_counts.items():
                        inclusion_n, exon_type, read_length, max_junction_n = (
                            coverage[sample_name]
                        )
                        psi = calculate_psi(
                            inclusion_n, exclusion_n, exon_length, read_length
                        )
                        fields = [
                            exon_type,
                            exon.chrom,
                            exon.start,
                            exon.end,
                            exon.strand,
                            inclusion_n,
                            max_junction_n,
                            junction.start,
                            junction.end,
                            exclusion_n,
                            sample_name,
                            psi,
                        ]
                        output.write("\t".join(map(str, fields)) + "\n")

                if exon_index % 1000 == 0:
                    logging.info("Processed %s exons", f"{exon_index:,}")

        logging.info("Calculating PSI for introns")
        with pysam.TabixFile(options.jxn_file) as junctions:
            for junction_index, junction_text in enumerate(junctions.fetch()):
                junction = JunctionRecord(junction_text)
                if junction_index % 1000 == 0:
                    logging.info("Processed %s junctions", f"{junction_index:,}")
                if junction.length < 20:
                    continue

                coverage = {
                    sample_name: get_coverage(
                        junction.contig,
                        junction.start,
                        junction.end,
                        junction.strand,
                        bam,
                        options.strandedness,
                    )
                    for sample_name, bam in bam_list
                }
                for sample_name, exclusion_n in junction.sample_counts.items():
                    (
                        inclusion_n,
                        _,
                        read_length,
                        max_junction_n,
                    ) = coverage[sample_name]
                    psi = calculate_psi(
                        inclusion_n, exclusion_n, junction.length, read_length
                    )
                    fields = [
                        "intron",
                        junction.contig,
                        junction.start,
                        junction.end,
                        junction.strand,
                        inclusion_n,
                        max_junction_n,
                        ".",
                        ".",
                        exclusion_n,
                        sample_name,
                        psi,
                    ]
                    output.write("\t".join(map(str, fields)) + "\n")

    for _, bam in bam_list:
        bam.close()


if __name__ == "__main__":
    main()
