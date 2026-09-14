## SWARM benchmarking scripts

## Read-level benchmark

Download the benchmarking dataset (IVT-m_data.tar.gz) from zenodo: https://zenodo.org/uploads/22123294

```
# extract tar archive 
tar -xzf IVT-m_data.tar.gz

# run read-level models on IVT data and evaluate predictions against sample modification labels
bash run_read-level_benchmark.sh

```

## Site-level benchmark

### First run read-level and site-level models on a transcriptome dataset
```
# read-level
singularity exec --nv tensorflow_24.01-tf2-py3-pysam.sif python3 SWARM_read_level.py -m $MOD --sam $SAM --fasta $FASTA --raw $BLOW5 -o $OUT

# cat to merge outputs if dataset has multiple replicates
cat rep1.pred.tsv rep2.pred.tsv > merged.pred.tsv

# sort read-level predictions 

sort -k 1 merged.pred.tsv > sorted.merged.pred.tsv

# site-level
singularity exec --nv tensorflow_24.01-tf2-py3-pysam.sif python3 SWARM_site_level.py -i sorted.merged.pred.tsv -o site-level.tsv
```
### Run liftover of transcriptomic to genomic coordinates
Install R2Dtool from : https://github.com/comprna/R2Dtool
```
# covert site-level output to bed
bash /PATH/TO/R2Dtool/scripts/cheui_to_bed.sh site-level.tsv site-level.bed

# run liftover (requires gtf matching the reference used for previous alignment steps)
r2d liftover -H -g human.gtf -i site-level.bed > site-level.lifted.bed
```
### Run site-level benchmark

#### Single-nucleotide coordinates (regular datasets)
```
get_validated_genomic_prc.py -v validated.bed -i site-level.lifted.bed -o output.prc.tsv
```

#### Range coordinates (pseU bisulfite datasets)
```
get_validated_genomic-range_prc.py validated.range.bed site-level.lifted.bed output.prc.tsv
```
