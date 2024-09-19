#!/usr/bin/python
import sys
import pandas as pd
import pysam 
import re 
import numpy as np
from collections import defaultdict, Counter
from scipy.stats import betabinom
import optparse
import logging



def beta_binomial_llr(grp1_lst, grp2_lst):

    mean_psi_grp1 = np.nanmean([x[3] for x in grp1_lst])
    mean_psi_grp2 = np.nanmean([x[3] for x in grp2_lst])

    stdv_psi_grp1 = np.nanstd([x[3] for x in grp1_lst])
    stdv_psi_grp2 = np.nanstd([x[3] for x in grp2_lst])

    alpha_grp1 = mean_psi_grp1/stdv_psi_grp1
    beta_grp1 = (1 - mean_psi_grp1)/stdv_psi_grp1

    alpha_grp2 = mean_psi_grp2/stdv_psi_grp2
    beta_grp2 = (1 - mean_psi_grp2)/stdv_psi_grp2

    llr_a = 0

    for (psi, inclusion_n, exclusion_n, mean) in grp1_lst:
        n = inclusion_n + exclusion_n
        p = betabinom.logpmf(inclusion_n, n, alpha_grp1, beta_grp1, loc = 0)
        llr_a += p

    for (psi, inclusion_n, exclusion_n, mean) in grp2_lst:
        n = inclusion_n + exclusion_n
        p = betabinom.logpmf(inclusion_n, n, alpha_grp2, beta_grp2, loc = 0)
        llr_a += p

    mean_psi_null = np.nanmean([x[3] for x in grp1_lst + grp2_lst])
    stdv_psi_null = np.nanstd( [x[3] for x in grp1_lst + grp2_lst])
    alpha_null = mean_psi_null/stdv_psi_null
    beta_null = (1 - mean_psi_null)/stdv_psi_null

    llr_null = 0
    for (psi, inclusion_n, exclusion_n, mean) in grp1_lst + grp2_lst:
        n = inclusion_n + exclusion_n
        p = betabinom.logpmf(inclusion_n, n, alpha_null, beta_null, loc = 0)
        llr_null += p

    return llr_a - llr_null, mean_psi_grp1, mean_psi_grp2, stdv_psi_grp1, stdv_psi_grp2


def main():

    args = sys.argv
    logging.basicConfig(level = logging.INFO, format = '%(asctime)s %(levelname)s:\t%(message)s')

    usage_instructions = f"""CryEx_psi_calculator -f fofn -e exon_coordinates -o output_file"""

    parser = optparse.OptionParser(usage = usage_instructions)

    parser.add_option('-f', dest='fofn', \
        help = 'Metadata file [required]', default = None)
    parser.add_option('-p', dest='input_PSI_file', \
        help = 'Input PSI file [required]')
    parser.add_option('-o', dest='output_file', \
        help = 'Output File')

    options, args = parser.parse_args()


    ### read metadata
    METADATA = pd.read_csv(options.fofn, sep = '\t')
    unique_group_values = METADATA['GROUP'].unique()

    assert len(unique_group_values) == 2, "Only two groups are allowed"

    GROUPS = {}
    for row in METADATA.itertuples():
        GROUPS[row.SAMPLE] = row.GROUP


    ### read psi file
    EXONS = defaultdict(lambda: defaultdict(list))
    psi = pd.read_csv(options.input_PSI_file, sep = '\t')
    # outf.write(f"exon_type\tchrom\texon_3ss\texon_5ss\tstrand\tinclusion_n\texc_5ss\texc_3ss\texclusion_n\tSAMPLE\tPSI\n")
    for row in psi.itertuples():
        exon  = (row.chrom, row.exon_3ss, row.exon_5ss, row.strand)
        group = GROUPS[row.SAMPLE]
        
        try: count_mean = row.inclusion_n/(row.inclusion_n + row.exclusion_n)
        except: count_mean = np.nan

        EXONS[exon][group].append((row.PSI, row.inclusion_n, row.exclusion_n, count_mean))


    ### diff splicing
    GRP1, GRP2 = unique_group_values
    for exon, group_dict in EXONS.items():
        llr, mean_psi1, mean_psi2, stdv_psi1, stdv_psi2 = beta_binomial_llr(group_dict[GRP1], group_dict[GRP2])
        psi_grp1_lst = ','.join(map(str, [round(x[0], 2) for x in group_dict[GRP1]]))
        psi_grp2_lst = ','.join(map(str, [round(x[0], 2) for x in group_dict[GRP2]]))

        lab  = f"{GRP1}-{GRP2}"
        exon = "{}:{}-{}:{}".format(*exon)
        print(exon, round(llr, 3), lab, round(mean_psi1 - mean_psi2, 2), psi_grp1_lst, psi_grp2_lst, sep = '\t')



        

    logging.info(f'Calculating psi of introns')

	



if __name__ == "__main__":
	main()
