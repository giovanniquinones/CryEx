# CryEx

A Python pipeline to identify Cryptix exons from RNA-Seq data.
It uses Stringtie to assemble transcripts and identify cryptic exons. 
It then calculates the PSI of each cryptic exons


## Installation

( we should replace this with the real package url)
python3 -m pip install --index-url https://test.pypi.org/simple/ --no-deps example-package-YOUR-USERNAME-HERE


### Dependencies
- stringtie
- multiprocess
- numpy
- pandas
- pysam
- subprocess

## Usage

```
# identify cryptic and annotated exons
Cryex_stringtie -f ${FOFN.tsv} -o ${EXONS.GTF}

# calculate splice junction usage
Cryex_junctions -f ${FOFN.tsv} -o ${JXN.BED}

# calculate PSI
Cryex_psi_calculator -f ${FOFN.tsv} -e ${EXONS.GTF} -j ${JXN.BED.GZ} -o {PSI.TSV} 
```

## Input

FOFN should be tab separated and have the following columns:
```
SAMPLE	BAM	STAR_SJ_OUT
sample1    /path/to/sample1.bam	/path/to/sample1.SJ.out.tab
sample2    /path/to/sample2.bam	/path/to/sample2.SJ.out.tab
sample3    /path/to/sample3.bam	/path/to/sample3.SJ.out.tab
```

## Output

```
exon_type       chrom   exon_3ss    exon_5ss    strand  inclusion_n     exc_5ss exc_3ss exclusion_n SAMPLE  PSI
first_exon      chr21   9907191     9907492     -       97              9896772 9966321 1           r2      0.96
first_exon      chr21   9907191     9907492     -       67              9896772 9966321 0           r3      1.0
first_exon      chr21   9907191     9907492     -       99              9896772 9966321 0           r4      1.0
first_exon      chr21   9907191     9907492     -       97              9896772 9966321 3           r2      0.92
first_exon      chr21   9907191     9907492     -       67              9896772 9966321 0           r3      1.0
first_exon      chr21   9907191     9907492     -       99              9896772 9966321 0           r4      1.0
```