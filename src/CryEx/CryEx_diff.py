#!/usr/bin/env python3
import logging
import optparse
import warnings
from collections import defaultdict
from collections.abc import Sequence

import numpy as np
import pandas as pd
from scipy.stats import betabinom, chi2
from statsmodels.stats.multitest import multipletests


SampleCounts = tuple[float, int, int, int]

np.seterr(all="ignore")
warnings.filterwarnings("ignore")


def mom_params(values: Sequence[float]) -> tuple[float, float, float]:
    """Estimate beta-distribution parameters using the method of moments."""
    mean = np.nanmean(values)
    variance = np.nanvar(values)
    phi = mean * (1 - mean) / variance - 1
    phi = max(phi, 1e-6)  # Guard against variance >= mean * (1 - mean).
    return variance, mean * phi, (1 - mean) * phi


def beta_binomial_llr(
    group1: list[SampleCounts], group2: list[SampleCounts]
) -> tuple[float, float, float]:
    """Return minimum group variance, log-likelihood ratio, and p-value."""
    psi1 = [sample[1] / sample[3] for sample in group1]
    psi2 = [sample[1] / sample[3] for sample in group2]

    variance1, alpha1, beta1 = mom_params(psi1)
    variance2, alpha2, beta2 = mom_params(psi2)

    alternative_log_likelihood = 0.0
    for _, inclusion_n, _, total_n in group1:
        alternative_log_likelihood += betabinom.logpmf(
            inclusion_n, total_n, alpha1, beta1
        )
    for _, inclusion_n, _, total_n in group2:
        alternative_log_likelihood += betabinom.logpmf(
            inclusion_n, total_n, alpha2, beta2
        )

    variance_null, alpha_null, beta_null = mom_params(psi1 + psi2)
    null_log_likelihood = 0.0
    for _, inclusion_n, _, total_n in group1 + group2:
        null_log_likelihood += betabinom.logpmf(
            inclusion_n, total_n, alpha_null, beta_null
        )

    log_likelihood_ratio = alternative_log_likelihood - null_log_likelihood
    p_value = 1 - chi2.cdf(2 * log_likelihood_ratio, 2)
    return (variance_null, log_likelihood_ratio, p_value)


def main() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(filename)s %(asctime)s %(levelname)s:\t%(message)s",
    )
    usage = "CryEx_diff -f fofn -p input_psi_file -o output_file"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-f", dest="fofn", help="Metadata file [required]")
    parser.add_option("-p", dest="input_PSI_file", help="Input PSI file [required]")
    parser.add_option("-o", dest="output_file", help="Output file [required]")
    parser.add_option(
        "--min_delta_psi",
        dest="min_delta_psi",
        help="Minimum absolute PSI difference [default=0.1]",
        default=0.1,
        type=float,
    )
    parser.add_option(
        "--min_llr",
        dest="min_llr",
        help="Minimum log-likelihood ratio [default=2.0]",
        default=2.0,
        type=float,
    )
    parser.add_option(
        "--min_cov",
        dest="min_cov",
        help="Minimum coverage per sample [default=10]",
        default=10,
        type=int,
    )
    parser.add_option(
        "--fdr",
        dest="fdr",
        help="False discovery rate threshold [default=0.05]",
        default=0.05,
        type=float,
    )
    options, _ = parser.parse_args()

    if not options.input_PSI_file or not options.fofn or not options.output_file:
        parser.error(f"Usage: {usage}")

    metadata = pd.read_csv(options.fofn, sep="\t")
    group_values = metadata["GROUP"].unique()
    if len(group_values) != 2:
        parser.error("Exactly two groups are required")

    sample_groups = dict(zip(metadata["SAMPLE"], metadata["GROUP"]))
    exon_counts: defaultdict[
        tuple[object, ...], defaultdict[object, list[SampleCounts]]
    ] = defaultdict(lambda: defaultdict(list))

    psi_dataframe = pd.read_csv(options.input_PSI_file, sep="\t")
    for row in psi_dataframe.itertuples():
        exon_key = (
            row.exon_type,
            row.chrom,
            row.exon_3ss,
            row.exon_5ss,
            row.strand,
            row.exc_5ss,
            row.exc_3ss,
        )
        group = sample_groups[row.SAMPLE]
        inclusion_n = row.jxn_inc_n + 1
        exclusion_n = row.exclusion_n + 1
        total_n = inclusion_n + exclusion_n
        exon_counts[exon_key][group].append(
            (row.PSI, inclusion_n, exclusion_n, total_n)
        )

    group1_name, group2_name = group_values
    logging.info("Processing %s regions", len(exon_counts))
    output_rows: list[list[object]] = []

    for exon_key, groups in exon_counts.items():
        group1 = groups[group1_name]
        group2 = groups[group2_name]
        if not group1 or not group2:
            continue

        sample_counts = group1 + group2
        if any(sample[3] < options.min_cov for sample in sample_counts):
            continue
        if any(np.isnan(sample[0]) for sample in sample_counts):
            continue

        variance, llr, p_value = beta_binomial_llr(group1, group2)
        if np.isnan(llr) or variance < 0.004:
            continue

        psi_group1 = ",".join(str(sample[0]) for sample in group1)
        psi_group2 = ",".join(str(sample[0]) for sample in group2)
        coverage_group1 = ",".join(str(sample[3]) for sample in group1)
        coverage_group2 = ",".join(str(sample[3]) for sample in group2)
        
        delta_psi = np.mean([sample[0] for sample in group1]) - np.mean(
            [sample[0] for sample in group2]
        )

        (
            exon_type,
            chrom,
            exon_start,
            exon_end,
            strand,
            junction_start,
            junction_end,
        ) = exon_key
        exon_coordinates = f"{chrom}:{exon_start}-{exon_end}:{strand}"
        flanking_junctions = f"{junction_start},{junction_end}"
        output_rows.append(
            [
                exon_type,
                exon_coordinates,
                flanking_junctions,
                round(llr, 3),
                p_value,
                np.nan,
                round(delta_psi, 2),
                f"{group1_name}-{group2_name}",
                psi_group1,
                psi_group2,
                coverage_group1,
                coverage_group2,
            ]
        )

    columns = [
        "exon_type",
        "exon_coords",
        "flanking_jxns",
        "LLR",
        "Pvalue",
        "Sig",
        "DeltaPSI",
        "DIFF",
        "PSIGroup1",
        "PSIGroup2",
        "CovGroup1",
        "CovGroup2",
    ]
    results = pd.DataFrame(output_rows, columns=columns)
    if not results.empty:
        _, adjusted_p_values, _, _ = multipletests(
            results["Pvalue"], method="fdr_bh"
        )
        results["Pvalue_Adj"] = adjusted_p_values
        results["Sig"] = (
            (results["LLR"] >= options.min_llr)
            & (results["DeltaPSI"].abs() >= options.min_delta_psi)
            & (results["Pvalue_Adj"] <= options.fdr)
        )
        results["Pvalue"] = results["Pvalue"].map(lambda value: f"{value:.2e}")
        results["Pvalue_Adj"] = results["Pvalue_Adj"].map(
            lambda value: f"{value:.2e}"
        )
    else:
        results["Pvalue_Adj"] = pd.Series(dtype="object")

    results.to_csv(options.output_file, sep="\t", index=False)
    logging.info("Differential splicing results written to %s", options.output_file)


if __name__ == "__main__":
    main()
