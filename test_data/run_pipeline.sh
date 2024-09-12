#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y # Error stream is merged with the standard output
#$ -l h_data=66G,h_rt=2:00:00
#$ -pe shared 1
#$ -r n # job is NOT rerunable
#$ -o SGE_STDOUT.txt 


# Input
FOFN=fofn.txt

# Output
EXONS=out/exons.gtf
JXN=out/jxn.bed
PSI=out/psi.txt

#================================================================================

# Identify cryptic and annotated exons
CryEx_stringtie -f ${FOFN} -o ${EXONS}

# Calculate splice junction usage
CryEx_junctions -f ${FOFN} -o ${JXN}

# Sort and gz compress bed file 
less ${JXN} | sort -k1,1V -k2,2n -k3,3n | bgzip > ${JXN}.gz

# Created tabix index
tabix -p bed -f ${JXN}.gz

# Calculate PSI
CryEx_psi_calculator -f ${FOFN} -e ${EXONS} -j ${JXN}.gz -o ${PSI} 
