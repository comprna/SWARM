# SWARM: Single-molecule Workflow for Analysing RNA Modifications 
Detection of pseudouridine, m6A, m5C, and ac4c on individual molecules from direct RNA nanopore signals


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
     * [Read-level detection](#read-level-detection)
     * [Site-level detection](#site-level-detection)
   * [Train new models](train-new-models)

------------------------------------------
# Dependencies
------------------------------------------
```
python=3.7
numpy==1.19.2
pandas==1.3.4
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
We highly recommend opting for f5c and sam files.
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
#build and compile htslib and slow5tools
bash build.sh

#compile SWARM executable
make
```

## Read-level detection

### sam + slow5 preprocessing (preferred)

Our preferred approach relies on gnu parallel which should be available on most HPC systems. 
Alternatively, preprocessing and prediction can be run separately, but that involves large temp files.


Example bash code to run SWARM read-level prediction.

```
module load tensorflow/2.8.0
module load parallel


predict_m1 () {
        SCRIPT=predict_model1_SWARM_bin_parallel.py
        MODEL=Model_100_epoch_relu.h5
        PATH_OUT=HepG2_mRNA_IVT_rep1.pred.tsv
        PATH_IN=temp.cpp.binh1
        python3 $SCRIPT -i $PATH_IN -o $PATH_OUT -m $MODEL -l Hep2G_1
}
export -f predict_m1


preprocess_cpp () {
        SCRIPT_CPP=SWARM2.15
        SAM=HepG2_mRNA_IVT_rep1.event.sam
        BLOW5=HepG2_IVT_rep1_fast5.blow5
        FASTA=Homo_sapiens.GRCh38.cdna.all.fa
        MODEL_KMER=IVT_model_c++.csv
        OUT=temp.cpp.binh1
        BASE=T
        $SCRIPT_CPP --sam $SAM --raw $BLOW5 --fasta $FASTA -m $MODEL_KMER -o $OUT --base $BASE
}

export -f preprocess_cpp


parallel ::: preprocess_cpp predict_m1
```

### eventalign.tsv preprocessing



## Site-level detection

## Mod-bam?

------------------------------------------
# Train new models
------------------------------------------
