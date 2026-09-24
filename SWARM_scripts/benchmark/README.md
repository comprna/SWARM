## SWARM benchmarking scripts

## Read-level benchmark
Benchmarking dataset and environment container are available on zenodo: https://zenodo.org/uploads/22123294

```
cd SWARM/SWARM_scripts/benchmark/read_level

# download the container and benchmarking dataset
wget https://zenodo.org/records/22123294/files/SWARM.sif
wget https://zenodo.org/records/22123294/files/IVT-m_data.tar.gz

# extract tar archive 
tar -xzf IVT-m_data.tar.gz

# run read-level models on IVT data (~1h total for all models and samples on 1 GPU Volta)
bash run_read-level_benchmark.sh

# should produce:
#  read-level .tsv files in outputs/
#  precision/recall curves (.png) for each modification combination

```

## Site-level benchmark
Raw reads obtained for this study are available from https://www.ebi.ac.uk/ena/browser/view/PRJEB106788 and public data from other studies is outlined in the paper

Instructions on processing raw reads for SWARM outputs are outlined at https://github.com/comprna/SWARM/tree/main#preprocess-raw-signals
### run read-level and site-level models on a transcriptome dataset
```
# read-level
singularity exec --nv SWARM.sif python3 /opt/SWARM/SWARM_scripts/SWARM_read_level.py -m $MOD --sam $SAM --fasta $FASTA --raw $BLOW5 -o $OUT

# cat to merge outputs if dataset has multiple replicates
cat rep1.pred.tsv rep2.pred.tsv > merged.pred.tsv

# sort read-level predictions 
sort -k 1 merged.pred.tsv > sorted.merged.pred.tsv

# site-level
singularity exec --nv SWARM.sif python3 /opt/SWARM/SWARM_scripts/SWARM_site_level.py -i sorted.merged.pred.tsv -o site-level.tsv
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
