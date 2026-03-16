# Preprocessing raw data

## Basecalling
SWARM was trained and tested with following basecaller versions and models, use the same versions for results comparable to our benchmarks:

| sequencing kit | basecaller version | basecaller model         |
|----------------|--------------------|--------------------------|
| sqk-RNA002     | guppy 6.4.6        | rna_r9.4.1_70bps_hac |
| sqk-RNA004     | dorado 0.7.2       | rna004_130bps_sup@v5.0.0 |

Recommended parameters RNA002:
```bash
guppy_basecaller -i $INPUTDIR --recursive -s $output_path.fastq -c guppy/ont-guppy/data/rna_r9.4.1_70bps_hac.cfg --device cuda:all:100%
```

Recommended parameters RNA004:
```bash
MODEL=dorado-0.7.2-linux-x64/rna004_130bps_sup@v5.0.0
dorado basecaller $MODEL $INPUTDIR -r -x cuda:all --emit-fastq > $output_path.fastq
```

## Alignment
minimap 2.24 for alignment and samtools 1.22 for quality control

Recommended parameters:
-k 5 for sythetic IVTs and -k 14 for human transcriptomes 

Example transcriptome alignment:

```bash
minimap2 -ax map-ont -k 14 ${fasta} ${input_path}/guppy_pass.fastq | samtools sort -o ${output_path}.bam
samtools index ${output_path}.bam

samtools view -b -F 2324  ${bam_file}.bam > ${bam_file}_pass_filtered.bam
samtools index ${bam_file}_pass_filtered.bam
```

## fast5 to slow5
This step is highly recommended, especially for large datasets.

Install slow5tools from:  https://github.com/hasindu2008/slow5tools

Example conversion command:

```bash
#convert fast5 files to slow5 files using 8 I/O processes
slow5tools f2s $INPUT_DIR -d $TEMPDIR  -p 8

#Merge all the slow5 files in to a single file using 8 threads
slow5tools merge $TEMPDIR -o $OUTDIR/${SAMPLE}.blow5 -t 8

#remove the temporary directory
rm -rf  $TEMPDIR
```

## Event alignment

### f5c
Our workflow supports both f5c .sam and nanopolish .tsv formats.
We highly recommend opting for **f5c** and **sam** files.
This requires the slow5 conversion outlined in previous step.

f5c is available from:  https://github.com/hasindu2008/f5c

Example event align command:

```bash
f5c index -t 48 $FASTQ_PATH --slow5 $SLOW5_PATH

f5c eventalign -t 48  -r $FASTQ_PATH --rna  -g $genome -b $BAM --slow5 $SLOW5_PATH --min-mapq 0 --signal-index --scale-events --samples --print-read-names --sam > $OUT
```

### nanopolish
We used this format in earlier stages of the project, our workflow can still support it.
Note that our prediction workflow is optimised for f5c sam format.

nanopolish f5c is available from:  https://github.com/jts/nanopolish

Example event align command:

```bash
nanopolish index -d ${fast5_path} -s ${guppy_files}/sequencing_summary.txt $fastq

nanopolish eventalign -t 48 --reads $fastq --bam $bam_file \
        --genome $fasta --signal-index --scale-events --samples --print-read-names > $output_path
```
