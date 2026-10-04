
# Input
set -euo pipefail

FOFN=fofn.txt

# Output
EXONS=out/exons.gtf
JXN=out/jxn.bed
PSI=out/psi.txt
DIFF=out/diff.txt

mkdir -p out

#================================================================================


# Identify cryptic and annotated exons
time CryEx_stringtie -f "${FOFN}" -o "${EXONS}"

# Calculate splice junction usage
time CryEx_junctions -f "${FOFN}" -o "${JXN}"

# Calculate PSI
time CryEx_psi_calculator -f "${FOFN}" -e "${EXONS}" -j "${JXN}.gz" -o "${PSI}"

# calculate diffential splicing
time CryEx_diff -f "${FOFN}" -p "${PSI}" -o "${DIFF}"
