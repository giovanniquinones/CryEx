#!/usr/bin/python
import sys
import pandas as pd
import os
import subprocess
import logging
import optparse


def main():
    args = sys.argv

    logging.basicConfig(level = logging.INFO, format = '%(asctime)s %(levelname)s:\t%(message)s')

    # option parser
    usage_instructions = f"""CryEx_stringtie -f fofn -o output_file"""

    parser = optparse.OptionParser(usage = usage_instructions)
    parser.add_option('-f', dest='fofn', \
        help = 'Metadata file [required]', default = None)
    parser.add_option('-s', dest='strandedness', \
        help = 'Library strandedness [optional]', default = 'fr-firststrand')
    parser.add_option('-o', dest='output', \
        help = 'Output File [required]', default = None)

    options, args = parser.parse_args()

    logging.info(f'bam file strandedness = {options.strandedness}')

    METADATA = pd.read_csv(options.fofn, sep = '\t')

    args = [(row.BAM, row.STAR_SJ_OUT, row.SAMPLE) for row in METADATA.itertuples()]

    bams = ' '.join(METADATA.BAM.values)

    STRINGTIE_COMMAND = f"stringtie {bams} -j 2 -f 0.05 -g 10 > {options.output}"

    mycmd = subprocess.getoutput(STRINGTIE_COMMAND)

    logging.info(mycmd)

if __name__ == "__main__":
    main()
