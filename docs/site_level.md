# Site-level detection

SWARM_site_level.py aggregates read-level data per reference site to compute:
 
 1. Stoichiometry (modification rate) computed from read-level modification calls
 
 2. Site-level probability (can narrow down modified sites from millions of tested coordinates)


## Sort read-level output

Read-level predictions **must be sorted** to ensure correct site-level aggregation.

```bash
# Use cat if pooling multiple replicates
cat rep1.pred.tsv rep2.pred.tsv > reads.pred.tsv

# Sort based on the first column
sort -k 1 reads.pred.tsv > reads.pred.tsv.sorted

# Remove temp file if pooling
rm -f reads.pred.tsv
```


## Running site-level detection

Run SWARM_site_level.py on sorted read-level data: 

```bash
python3 SWARM_site_level.py -i <pred.tsv.sorted> -o <OUT> -d <0.5,0.5>
```

```bash
Required:
  -i, --input           Path to the sorted read-level prediction file
  -o, --file_out        Path to the output file

Optional:
  -d, --double_cutoff   Read-level cutoffs for computing stoichiometry [0.5,0.5]
  -n, --min_reads       Minimun read coverage to output a site [20]
  -c, --cutoff          Site-level probability cutoff for printing sites [0.0]
  -m, --DL_model        Custom path to pretrainned site-level model
  --arch                NN architecture [Mini/Mid/Large]
  -h, --help            Show this help message and exit

```


## Stoichiometry parameters

` -d, --double_cutoff ` arg can be used to select read-level cutoffs for computing stoichiometry. 

Expected input is a ',' separated value providing cutoffs for unmodified and modified calls. 

Default value is 0.5,0.5 ; meaning that p < 0.5 is unmodified and p > 0.5 is modified. 

Values between the two cutoffs are not included in stoichiometry computation, but are still kept for site-level model prediction and included in the coverage. 

Using higher cutoffs for calling modified reads would reduce stoichiometry but should give lower false positive rates.

Stoichiometry is reported as rate (between 0 and 1).


## Output format

SWARM_site_level.py outputs a tsv file in the same fromat as CHEUI (Mateos et al., 2024):
```bash
contig	position	site	coverage	stoichiometry	probability
ENST00000000233.10	1000	GATCTTGAG	223	0.04522613065326633	4.0531158447265625e-06
ENST00000000233.10	1001	ATCTTGAGT	227	0.06111111111111111	8.106231689453125e-06
ENST00000000233.10	1012	TAAATTTGC	164	0.08527131782945736	6.079673767089844e-06
ENST00000000233.10	1013	AAATTTGCT	182	0.07586206896551724	4.1484832763671875e-05
ENST00000000233.10	1017	TTGCTGTGG	151	0.00684931506849315	0.00010669231414794922
```

## Post-processing 

### Site-level cutoffs

SWARM can detect modifications in a single sample, where the target modfiication is usually present in under 1% of the tested coordinates. 

To enrich for the true modified sites, we use a standard 10% stoichiometry cutoff with additional neural-network filtering based on the distribution of read-level probabilities.  

Site-level models were trained on pre-defined mixtures of modified and unmodified reads and extensively benchmarked on IVT transcriptomes and cellular data. 

For single-sample detection we select sites with: 

1. Stoichiometry > 0.1 (>10%)

2. Site-level probability with a false-positive rate < 0.1% on unmodified IVT transcriptomes

Cutoffs for different contexts:

| Kit    | Modification | Context         | Cutoff          |
|--------|--------------|-----------------|-----------------|
| RNA002 | m6A          | all-context     | 0.999953        |
| RNA002 | m6A          | DRACH-only      | 0.9972          |
| RNA004 | m6A          | all-context     | 0.9999999999986 |
| RNA004 | m6A          | DRACH-only      | 0.9999986       |
| RNA002 | Ψ            | all-context     | 0.9943          |
| RNA002 | Ψ            | PUS7-TRUB1-only | 0.9925          |
| RNA004 | Ψ            | all-context     | 0.999999988     |
| RNA004 | Ψ            | PUS7-TRUB1-only | 0.999999975     |
| RNA002 | m5C          | all-context     | 0.999969        |
| RNA002 | m5C          | NSUN6-only      | 0.99981         |
| RNA004 | m5C          | all-context     | 0.99999999989   |
| RNA004 | m5C          | NSUN6-only      | 0.988           |

### Motif filtering

### Liftover to genome
