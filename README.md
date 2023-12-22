# SWARM: Single-molecule Workflow for Analysing RNA Modifications 
Detection of pseudouridine, m6A, m5C, and ac4C on individual molecules from direct RNA nanopore signals


## About


------------------------------------------
# Table of Contents
------------------------------------------

   * [Dependencies](#dependencies)
   * [Preprocess raw signals](#preprocess-raw-signals)
     * [Basecalling](#basecalling)
     *   [Alignment](#alignment)
     *   [fast5 to slow5](#fast5-to-slow5)
     *   [Event alignment](#event-alignment)
   
   * [Detect RNA modificatios](#detect-rna-modifications)
     * [Installation](#installation)
     * [Read-level single-base detection](#read-level-single-base-detection)
     * [Site-level detection](#site-level-detection)
   * [Train new models](train-new-models)

------------------------------------------
# Dependencies
------------------------------------------
```
python=3.9
numpy==1.19.2
tensorflow-gpu==2.4.1
keras-preprocessing==1.1.2
```


------------------------------------------
# Preprocess raw signals
------------------------------------------

## Basecalling
We trained on signals basecalled with guppy 6.4.6 


Recommended parameters:
```
guppy_basecaller -i $INPUTDIR --recursive -s $output_path -c guppy/ont-guppy/data/rna_r9.4.1_70bps_hac.cfg --device cuda:all:100% --compress_fastq --gpu_runners_per_device 2
```


## Alignment
We used minimap 2.24 for alignment and samtools 1.12 for quality checks

Recommended parameters:
-k 5 for curlcakes and -k 14 for standard human transcriptomes 

```
minimap2 -ax map-ont -k 5 ${fasta} ${input_path}/guppy_pass.fastq | samtools sort -o ${output_path}.bam
samtools index ${output_path}.bam

samtools view -b -F 2324  ${bam_file}.bam > ${bam_file}_pass_filtered.bam
samtools index ${bam_file}_pass_filtered.bam
```


## fast5 to slow5
This step is optional but highly recommended, especially for large datasets.

https://github.com/hasindu2008/slow5tools

Example conversion command:

```
#convert fast5 files to slow5 files using 8 I/O processes
slow5tools f2s $INPUT_DIR -d $TEMPDIR  -p 8

#Merge all the slow5 files in to a single file using 8 threads
slow5tools merge $TEMPDIR -o $OUTDIR/${SAMPLE}.blow5 -t 8

#remove the temporary directory
rm -rf  $TEMPDIR
```


## Event alignment

### f5c

Our workflow supports both f5c sam and nanopolish tsv formats.
We highly recommend opting for **f5c** and **sam** files.
This requires the slow5 conversion outlined in previous step.

https://github.com/hasindu2008/f5c

Example event align command:

```
#f5c index
f5c index -t 48 $FASTQ_PATH --slow5 $SLOW5_PATH

#sam event alignment format
f5c  eventalign -t 48  -r $FASTQ_PATH --rna  -g $genome -b $BAM --slow5 $SLOW5_PATH --min-mapq 0 --secondary=yes --signal-index --scale-events --samples --print-read-names --sam > $OUT
```

### nanopolish

We used this format in earlier stages of the project, our workflow can still support it.
Note that our workflow is optimised for f5c sam format.

https://github.com/jts/nanopolish

Example event align command:

```
nanopolish index -d ${fast5_path} -s ${guppy_files}/sequencing_summary.txt $fastq

nanopolish eventalign -t 48 --reads $fastq --bam $bam_file \
        --genome $fasta --signal-index --scale-events --samples --print-read-names > $output_path
```



------------------------------------------
# Detect RNA modifications
------------------------------------------

## Installation


This step is highly recommended and required for using our C++ preprocessing of sam event files. Skip if using eventalign.tsv format.


```
cd SWARM_scripts/preprocessing/

#build and compile htslib, slow5tools, SWARM_preprocess
bash build.sh

```

## Read-level single-base detection

### sam + slow5 preprocessing (preferred)

Use this approach for faster and simultaneous preprocessing + model inference. Run build.sh from above section.

Tensorflow (GPU-configured) should be available on most HPC systems with GPU access. 
Otherwise, it is highly advised to use tensorflow configured for GPU. https://www.tensorflow.org/install/

Example bash code to run SWARM read-level prediction.

```
module load tensorflow


export MOD=pU
export FASTA=Homo_sapiens.GRCh38.cdna.fa
export RAW=Hek293_mRNA.blow5
export SAM=Hek293_mRNA_f5C.sam
export OUT=Hek293_mRNA_pU
export TEMP=$PBS_JOBFS     #for optional speed up provide SSD path if available ( such as jobfs on PBS systems)

python3 ./SWARM_scripts/SWARM_read_level.py -m $MOD --sam $SAM --fasta $FASTA --raw $RAW -o $OUT --temp $TEMP
```

### eventalign.tsv preprocessing

Alternatively, preprocessing and prediction can be run separately from eventalign.tsv, but that involves large temp files.

First preprocess the event alignments.

```

export MOD=pU
export BAM=Hek293_mRNA_f5C.bam
export EVENTS=Hek293_mRNA.events.tsv
export OUT=Hek293_mRNA_pU

python3 ./SWARM_scripts/SWARM_read_level.py --preprocess -m $MOD --bam BAM --nanopolish $EVENTS -o $OUT
```

Then predict modification states.

```

export MOD=pU
export PICKLE=Hek293_mRNA_pU_T.pickle
export OUT=Hek293_mRNA_pU.pred.tsv

python3 ./SWARM_scripts/SWARM_read_level.py --predict -m $MOD --pickle $PICKLE -o $OUT
```


## Site-level detection

## Mod-bam?

------------------------------------------
# Train new models
------------------------------------------
## Train read-level prediction
### Trim eventalign files
This optional step reduces the time to retrain models as preprocessing only a fraction of signals from a whole sample is usually enough for training. We trim for events comprising 500 signals per 9mer. 

```
python3 ./SWARM_scripts/train_models/trim_tsv_events.py -i <eventalign.tsv> -o <out_prefix> --limit-out 500
```
### Preprocess trimmed files
Preprocess trimmed files for model1 input features. Make sure to include **--out_counter** arg here!
```
python3 ./SWARM_scripts/SWARM_read_level.py --preprocess -m <pU/m6A/m5C/ac4C> --bam <BAM> //
--nanopolish <eventalign_trimmed.tsv> -o <out_prefix> --out_counter
```
### Split training/validation/testing data
 Use this step for stratified sampling of the preprocessed signals. Splits equal number of signals per 9mer for positive/negative labels in each train/test/validation set (60/20/20 split by default).
```
python3 ./SWARM_scripts/train_models/split_training_by_9mers.py -i <preprocesed.pickle> //
--counts <preprocessed.counts> -o <outpath> --limit <signals_per_9mer>
```


