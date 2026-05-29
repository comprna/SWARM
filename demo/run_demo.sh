#!/bin/bash

module load tensorflow/2.15.0
module load samtools

FASTA=reference.fa
BLOW5=reads.blow5
SAM=events.sam

SCRIPTDIR=../SWARM_scripts

for MOD in "m6A" "pU"; do
	OUT=$PWD/$MOD
	python3 $SCRIPTDIR/SWARM_read_level.py -m $MOD --sam $SAM --fasta $FASTA --raw $BLOW5 -o $OUT --modsam
	
	sort -k1 $OUT.pred.tsv > $OUT.pred.tsv.sorted
	python3 $SCRIPTDIR/SWARM_site_level.py -i $OUT.pred.tsv.sorted -o $OUT.sites.pred.tsv

	# get significant sites
	python3 $SCRIPTDIR/postprocess/get_significant_sites.py -i $OUT.sites.pred.tsv -o $OUT.sig.sites.bed --kit RNA002 -m $MOD --context writer
	
	# filter modbam for significant sites
	samtools view -b $OUT.mod.sam > $OUT.mod.bam
	python3 $SCRIPTDIR/postprocess/filter_modbam_sites.py $OUT.mod.bam $OUT.sig.mod.bam $OUT.sig.sites.bed
		
done

# merge m6A and pU modbams, filter read-level probabilities (p>0.9) for cleaner visualisation
echo "pU.sig.mod.bam,pU" > merge.input.csv
echo "m6A.sig.mod.bam,m6A" >> merge.input.csv
python3 $SCRIPTDIR/postprocess/merge_modbams.py --merge merge.input.csv -o merged.sig.mod.bam -t 0.9

