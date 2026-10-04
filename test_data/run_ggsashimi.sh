
## ggsashimi plot

cat fofn.txt | grep -v SAMPLE | cut -f1,2,4 > input_bams.tsv

printf "red\norange" > palette.txt

GTF=out/exons.gtf

./ggsashimi.py \
	-b input_bams.tsv \
	-g $GTF \
	-c "chr21:45457573-45472322" \
	-M 3 -C 1 -O 1 -F pdf \
	--base-size=14 --height=1.5 --ann-height=3 --width=7 --fix-y-scale -P palette.txt

## Output
echo "Check output sashimi.pdf"  
