# CryEx

(This repository is private for now)

A Python pipeline to identify Cryptix exons from RNA-Seq data.
It uses Stringtie to assemble transcripts and identify cryptic exons. 
It then calculates the PSI of each cryptic exons

Instead of processing one bam file at the time, it takes as input a metadata
file (FOFN) and can processes multiple bam files at once. 


## Installation
Download this repository then:

```
pip install dis/CryEx-0.0.1-py3-none-any.whl
```

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

# sort and gz compress bed file 
less ${JXN.BED} | sort -k1,1V -k2,2n -k3,3n | bgzip -c > ${JXN.BED.GZ}

# Created tabix index
tabux -p bed -f ${JXN.BED.GZ}

# calculate PSI
Cryex_psi_calculator -f ${FOFN.tsv} -e ${EXONS.GTF} -j ${JXN.BED.GZ} -o {PSI.TSV} 

# calculate diffential splicing
Cryex_diff -f ${FOFN.tsv} -p {PSI.TSV} -o {DIFF.tsv}
```

## Input

FOFN should be tab separated and have the following columns:
For differential splicing, Cryex will use the 'GROUP' column
```
SAMPLE	BAM	STAR_SJ_OUT	GROUP
sample1    /path/to/sample1.bam	/path/to/sample1.SJ.out.tab	KD
sample2    /path/to/sample2.bam	/path/to/sample2.SJ.out.tab	KD
sample3    /path/to/sample3.bam	/path/to/sample3.SJ.out.tab	CTRL
sample4	   /path/to/sample4.bam /path/to/sample4.SJ.out.tab     CTRL
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



