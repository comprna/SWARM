# Read-level modification calling

Use SWARM_read_level.py to detect modifications in individual molecules. 



## sam + slow5 preprocessing (default)

Use this approach for faster and simultaneous preprocessing + model inference. Must run build.sh from install section to compile C++ code for signal preprocessing.

Models for RNA002 or RNA004 chemistry are automatically selected by default based on the blow5 header. 

Example bash code to run SWARM read-level prediction for m6A:

```bash
module load tensorflow
conda activate SWARM
python3 SWARM_read_level.py -m <RNAmod> -s <SAM> -f <FASTA> -r <BLOW5> -o <OUT> 
```

```
 Required:
  -m RNAMOD, --RNAmod RNAMOD    Target RNA modification [m6A/pU/m5C]
  -o OUT, --out OUT             Path for the output tsv file
  -s SAM, --sam SAM             Path to the input sam event align
  -f FASTA, --fasta FASTA       Path to the input fasta reference genome
  -r RAW, --raw RAW             Path to the input signals in blow5 format

Optional:
  --model1 MODEL1               Custom path to the trained model1
  --kmer KMER                   Custom path to the kmer model
  --cpp CPP                     Custom path to the compiled c++ preprcessing binary
  --kit KIT                     RNA sequencing kit [RNA004/RNA002]
  --temp TEMP                   Directory for temp files
  --arch ARCH                   Model1 network, Mini is default. [Mini/Mid/Large].
  --modsam                      Outputs OUT.mod.sam (MM/ML tags) and OUT.pred.tsv
  --nworkers NWORKERS           Number of preprocessing workers (default 4, 1 for --modsam)
  -h, --help                    Show this help message and exit
```


## eventalign.tsv preprocessing

Alternatively, preprocessing and prediction can be run separately from eventalign.tsv, but that involves massive temp files (can be terabytes).

First preprocess the event alignments:

```bash
python3 SWARM_read_level.py --preprocess -m m6A --bam BAM --nanopolish $EVENTS -o $OUT.pickle
```

Then predict modifications:

```bash
python3 SWARM_read_level.py --predict -m m6A --pickle $OUT.pickle -o $OUT.pred.tsv
```


## Output format

### pred.tsv

Default output is in tsv format and contains 3 columns:

```bash
ENST00000390289.2_53_TCCCTCTCC_3f1d97d4-036e-478e-807a-dc68148e832d_326_28_T    0.01007  1
ENST00000390289.2_55_CCTCTCCCA_3f1d97d4-036e-478e-807a-dc68148e832d_328_37_T    0.03012  1
ENST00000390289.2_63_AGCCTGTGC_3f1d97d4-036e-478e-807a-dc68148e832d_336_42_T    0.00112  1
ENST00000390289.2_65_CCTGTGCTG_3f1d97d4-036e-478e-807a-dc68148e832d_338_42_T    0.00612  1
ENST00000390289.2_68_GTGCTGACT_3f1d97d4-036e-478e-807a-dc68148e832d_341_41_T    0.01945  1
```


1. Base metadata: refContig_refPosition_9mer_readID_readPosition_qscore_calledBase   

2. Molecule-level probability, between 0 and 1. Target modification should have higher probability.

3. Model code, used for selecting parameters for site-level prediction. 

```bash
model_code = {
    "pU_RNA002":    "1",
    "pU_RNA004":    "2",
    "m5C_RNA002":   "3",
    "m5C_RNA004":   "4",
    "m6A_RNA002":   "5",
    "m6A_RNA004":   "6"
}
```


### mod.sam

SWARM predictions can be encoded in modsam format using the --modsam tag 

Described in detail in **modification visualization** section. 